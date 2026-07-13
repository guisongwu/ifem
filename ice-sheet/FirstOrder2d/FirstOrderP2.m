function [soln,eqn,info,node,elem,bdFlag] = FirstOrderP2(option)
%% FIRSTORDERP2 two-dimensional x-z FO cross-section solver.
%
% Solve the y-independent FO/BP cross-section model on a slab
%
%   0 <= x <= L,  -slope*x <= z <= H-slope*x,
%
% where the second coordinate is z.  The model assumes v = 0 and solves for
% the horizontal velocity u(x,z) with quadratic triangular elements:
%
%   -d_x(4 eta u_x) - d_z(eta u_z) = f.
%
% Boundary conditions:
%   top z = H-slope*x: zero traction, natural boundary condition;
%   left/right x = 0,L: periodic, with vertical shift from the slope;
%   bottom z = -slope*x: linear sliding beta u, with m = 1 by default.

if nargin < 1, option = struct; end
option = setoption(option,'L',4);
option = setoption(option,'H',1);
option = setoption(option,'slope',0.1);
option = setoption(option,'h',[1,0.5]);
option = setoption(option,'tol',1e-8);
option = setoption(option,'maxIt',100);
option = setoption(option,'damping',0.8);
option = setoption(option,'eps_reg',1e-3);
option = setoption(option,'quadorder',4);
option = setoption(option,'edgequadorder',3);
option = setoption(option,'printlevel',1);
option = setoption(option,'periodic_x',[0,option.L]);
option = setoption(option,'periodic_tol',1e-10);
option = setoption(option,'assemble_tangent',false);
option = setoption(option,'residual_tol',option.tol);
option = setoption(option,'residual_check_threshold',...
    max(1e-2,sqrt(option.residual_tol)));

if isfield(option,'node') && isfield(option,'elem')
    node = option.node;
    elem = option.elem;
    if isfield(option,'bdFlag')
        bdFlag = option.bdFlag;
    else
        bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
    end
else
    box = [0,option.L,0,option.H];
    [node,elem] = rectmesh2(box,option.h);
    bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
    node(:,2) = node(:,2)-option.slope*node(:,1);
end

pde.A = getoption(option,'A',1);
pde.n = getoption(option,'n',3);
pde.rho = getoption(option,'rho',910);
pde.gravity = getoption(option,'gravity',9.81);
pde.beta = getoption(option,'beta',1e3);
pde.m = getoption(option,'m',1);
pde.gradS = getoption(option,'gradS',...
    @(p) -option.slope*ones(size(p,1),1));
pde.f = getoption(option,'f',[]);

[elem2dof,edge] = dofP2(elem);
[Dlambda,area] = gradbasis(node,elem);
N = size(node,1);
NT = size(elem,1);
NE = size(edge,1);
Ndof = N+NE;
dofNode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
[periodicRep,periodicInfo] = periodicrepresentatives(...
    dofNode,option.periodic_x,option.slope,option.periodic_tol);
[periodicMasterDof,~,scalarMaster] = unique(periodicRep);
Nmaster = length(periodicMasterDof);
Pperiodic = sparse((1:Ndof)',scalarMaster,1,Ndof,Nmaster);

masterState = zeros(Nmaster,1);
if isfield(option,'u0') && numel(option.u0) == Ndof
    masterState = accumarray(scalarMaster,option.u0(:),...
        [Nmaster,1],@mean);
end
state = Pperiodic*masterState;
residual = zeros(option.maxIt,1);
equationResidual = NaN(option.maxIt,1);
residualChecked = false(option.maxIt,1);
viscosityRange = zeros(option.maxIt,2);
checkResidualEveryStep = false;
t0 = cputime;

for k = 1:option.maxIt
    [A,etaMin,etaMax] = assembleviscous(state);
    [Ab,bedCoefficient] = assemblebed(state);
    F = assembleforce;
    A = A+Ab;
    Ared = Pperiodic'*A*Pperiodic;
    Fred = Pperiodic'*F;

    newMasterState = Ared\Fred;
    updatedMasterState = (1-option.damping)*masterState + ...
        option.damping*newMasterState;
    updatedState = Pperiodic*updatedMasterState;

    residual(k) = norm(updatedState-state)/max(1,norm(updatedState));
    masterState = updatedMasterState;
    state = updatedState;
    viscosityRange(k,:) = [etaMin,etaMax];

    if residual(k) <= option.residual_check_threshold
        checkResidualEveryStep = true;
    end
    if checkResidualEveryStep || k == option.maxIt
        [equationResidual(k),A,bedCoefficient,etaMin,etaMax] = ...
            evaluateresidual(state);
        viscosityRange(k,:) = [etaMin,etaMax];
        residualChecked(k) = true;
    end

    if option.printlevel > 0
        if residualChecked(k)
            fprintf(['FO 2D P2 %2d: relchange %.3e, residual %.3e, ',...
                'eta [%.3e, %.3e]\n'],...
                k,residual(k),equationResidual(k),etaMin,etaMax);
        else
            fprintf('FO 2D P2 %2d: relchange %.3e, eta [%.3e, %.3e]\n',...
                k,residual(k),etaMin,etaMax);
        end
    end
    if residualChecked(k) && residual(k) < option.tol && ...
            equationResidual(k) < option.residual_tol
        break
    end
end

[wDiagnostic,pDiagnostic] = diagnosticwp(state);
soln.u = state;
soln.w = wDiagnostic;
soln.p = pDiagnostic;
eqn.A = A;
eqn.edge = edge;
eqn.elem2dof = elem2dof;
eqn.dofNode = dofNode;
eqn.periodicProjection = Pperiodic;
eqn.periodicMasterDof = periodicMasterDof;

if option.assemble_tangent
    Kt = assembleviscoustangent(state);
    Kbt = assemblebedtangent(state);
    eqn.tangent = Pperiodic'*(Kt+Kbt)*Pperiodic;
    eqn.applyBetaDerivative = @assemblebetadirection;
end

info.itStep = k;
info.relchange = residual(1:k);
info.equationResidual = equationResidual(1:k);
info.residualChecked = residualChecked(1:k);
info.residualTolerance = option.residual_tol;
info.viscosityRange = viscosityRange(1:k,:);
info.converged = residualChecked(k) && residual(k) < option.tol && ...
    equationResidual(k) < option.residual_tol;
info.solveTime = cputime-t0;
info.bedCoefficient = bedCoefficient;
info.periodic = periodicInfo;
info.geometry.L = option.L;
info.geometry.H = option.H;
info.geometry.slope = option.slope;

if nargout == 0
    fprintf('FO 2D rectangular periodic solve: %d dofs, %d iterations\n',...
        length(soln.u),info.itStep);
    clear soln eqn info node elem bdFlag
end

    function [A,etaMin,etaMax] = assembleviscous(uk)
        [lambda,w] = quadpts(option.quadorder);
        rows = zeros(size(lambda,1)*36*NT,1);
        cols = zeros(size(lambda,1)*36*NT,1);
        vals = zeros(size(lambda,1)*36*NT,1);
        cursor = 0;
        etaMin = inf;
        etaMax = 0;
        for q = 1:size(lambda,1)
            Dphi = p2gradient(lambda(q,:),Dlambda);
            gradU = evalgrad(uk,elem2dof,Dphi);
            ux = gradU(:,1);
            uz = gradU(:,2);
            epsII = ux.^2+0.25*uz.^2;
            xq = evalpoint(lambda(q,:));
            Aq = coefficient(pde.A,xq);
            nq = coefficient(pde.n,xq);
            eta = 0.5.*Aq.^(-1./nq).*...
                (epsII+option.eps_reg^2).^((1-nq)./(2*nq));
            etaMin = min(etaMin,min(eta));
            etaMax = max(etaMax,max(eta));
            for a = 1:6
                da = Dphi(:,:,a);
                for b = 1:6
                    db = Dphi(:,:,b);
                    kab = w(q)*area.*eta.*...
                        (4*da(:,1).*db(:,1)+da(:,2).*db(:,2));
                    idx = cursor+(1:NT);
                    rows(idx) = elem2dof(:,a);
                    cols(idx) = elem2dof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        A = sparse(rows,cols,vals,Ndof,Ndof);
    end

    function [Ab,coefAtEdge] = assemblebed(uk)
        Ab = sparse(Ndof,Ndof);
        coefAtEdge = [];
        if ~any(bdFlag(:)==3), return; end
        [s,w] = gaussedge(option.edgequadorder);
        ii = [];
        jj = [];
        ss = [];
        coefAtEdge = zeros(nnz(bdFlag(:)==3),1);
        coefCounter = 0;
        for e = 1:3
            elemIdx = find(bdFlag(:,e)==3);
            if isempty(elemIdx), continue; end
            edgeDof = elem2dof(elemIdx,edgelocaldof(e));
            edgeLength = localedgelength(elemIdx,e);
            tangentX = localedgetangentx(elemIdx,e);
            for q = 1:length(w)
                lambdaEdge = edgelambda(e,s(q));
                phi = p2basis(lambdaEdge);
                edgePhi = phi(edgelocaldof(e));
                xq = evalpoint(lambdaEdge,elemIdx);
                beta = coefficient(pde.beta,xq);
                uq = uk(edgeDof)*edgePhi';
                tangentialU = uq.*tangentX;
                gamma = beta.*...
                    (tangentialU.^2+option.eps_reg^2).^...
                    ((pde.m-1)/2).*tangentX.^2;
                if q == ceil(length(w)/2)
                    coefAtEdge(coefCounter+(1:length(elemIdx))) = gamma;
                    coefCounter = coefCounter+length(elemIdx);
                end
                for a = 1:3
                    for b = 1:3
                        sAb = w(q)*edgeLength.*gamma*edgePhi(a)*edgePhi(b);
                        ii = [ii;edgeDof(:,a)]; %#ok<AGROW>
                        jj = [jj;edgeDof(:,b)]; %#ok<AGROW>
                        ss = [ss;sAb]; %#ok<AGROW>
                    end
                end
            end
        end
        Ab = sparse(double(ii),double(jj),ss,Ndof,Ndof);
    end

    function At = assembleviscoustangent(uk)
        [lambda,w] = quadpts(option.quadorder);
        rows = zeros(size(lambda,1)*36*NT,1);
        cols = zeros(size(lambda,1)*36*NT,1);
        vals = zeros(size(lambda,1)*36*NT,1);
        cursor = 0;
        for q = 1:size(lambda,1)
            Dphi = p2gradient(lambda(q,:),Dlambda);
            gradU = evalgrad(uk,elem2dof,Dphi);
            ux = gradU(:,1);
            uz = gradU(:,2);
            epsII = ux.^2+0.25*uz.^2;
            xq = evalpoint(lambda(q,:));
            Aq = coefficient(pde.A,xq);
            nq = coefficient(pde.n,xq);
            strain = epsII+option.eps_reg^2;
            exponent = (1-nq)./(2*nq);
            eta = 0.5.*Aq.^(-1./nq).*strain.^exponent;
            for a = 1:6
                da = Dphi(:,:,a);
                stateTest = 4*ux.*da(:,1)+uz.*da(:,2);
                for b = 1:6
                    db = Dphi(:,:,b);
                    strainProduct = 4*da(:,1).*db(:,1)+...
                        da(:,2).*db(:,2);
                    stateDirection = 2*ux.*db(:,1)+0.5*uz.*db(:,2);
                    kab = w(q)*area.*eta.*...
                        (strainProduct+...
                         exponent./strain.*stateTest.*stateDirection);
                    idx = cursor+(1:NT);
                    rows(idx) = elem2dof(:,a);
                    cols(idx) = elem2dof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        At = sparse(rows,cols,vals,Ndof,Ndof);
    end

    function Abt = assemblebedtangent(uk)
        Abt = sparse(Ndof,Ndof);
        if ~any(bdFlag(:)==3), return; end
        [s,w] = gaussedge(option.edgequadorder);
        ii = [];
        jj = [];
        ss = [];
        exponent = (pde.m-1)/2;
        for e = 1:3
            elemIdx = find(bdFlag(:,e)==3);
            if isempty(elemIdx), continue; end
            edgeDof = elem2dof(elemIdx,edgelocaldof(e));
            edgeLength = localedgelength(elemIdx,e);
            tangentX = localedgetangentx(elemIdx,e);
            for q = 1:length(w)
                lambdaEdge = edgelambda(e,s(q));
                phi = p2basis(lambdaEdge);
                edgePhi = phi(edgelocaldof(e));
                xq = evalpoint(lambdaEdge,elemIdx);
                beta = coefficient(pde.beta,xq);
                uq = uk(edgeDof)*edgePhi';
                tangentialU = uq.*tangentX;
                speed = tangentialU.^2+option.eps_reg^2;
                tangentCoefficient = beta.*...
                    tangentX.^2.*...
                    (speed.^exponent+(pde.m-1)*tangentialU.^2.*...
                    speed.^(exponent-1));
                for a = 1:3
                    for b = 1:3
                        sAb = w(q)*edgeLength.*tangentCoefficient*...
                            edgePhi(a)*edgePhi(b);
                        ii = [ii;edgeDof(:,a)]; %#ok<AGROW>
                        jj = [jj;edgeDof(:,b)]; %#ok<AGROW>
                        ss = [ss;sAb]; %#ok<AGROW>
                    end
                end
            end
        end
        Abt = sparse(double(ii),double(jj),ss,Ndof,Ndof);
    end

    function load = assemblebetadirection(betaDirection)
        fullLoad = zeros(Ndof,1);
        if ~any(bdFlag(:)==3)
            load = Pperiodic'*fullLoad;
            return;
        end
        [s,w] = gaussedge(option.edgequadorder);
        exponent = (pde.m-1)/2;
        for e = 1:3
            elemIdx = find(bdFlag(:,e)==3);
            if isempty(elemIdx), continue; end
            edgeDof = elem2dof(elemIdx,edgelocaldof(e));
            edgeLength = localedgelength(elemIdx,e);
            tangentX = localedgetangentx(elemIdx,e);
            for q = 1:length(w)
                lambdaEdge = edgelambda(e,s(q));
                phi = p2basis(lambdaEdge);
                edgePhi = phi(edgelocaldof(e));
                xq = evalpoint(lambdaEdge,elemIdx);
                deltaBeta = coefficient(betaDirection,xq);
                uq = state(edgeDof)*edgePhi';
                tangentialU = uq.*tangentX;
                tractionDirection = deltaBeta.*...
                    (tangentialU.^2+option.eps_reg^2).^exponent.*...
                    uq.*tangentX.^2;
                for a = 1:3
                    contribution = w(q)*edgeLength.*...
                        tractionDirection*edgePhi(a);
                    fullLoad = fullLoad+accumarray(...
                        double(edgeDof(:,a)),contribution,[Ndof,1]);
                end
            end
        end
        load = Pperiodic'*fullLoad;
    end

    function F = assembleforce
        F = zeros(Ndof,1);
        [lambda,w] = quadpts(option.quadorder);
        for q = 1:size(lambda,1)
            phi = p2basis(lambda(q,:));
            xq = evalpoint(lambda(q,:));
            if ~isempty(pde.f)
                fq = coefficient(pde.f,xq);
            else
                gradS = coefficient(pde.gradS,xq);
                rho = coefficient(pde.rho,xq);
                gravity = coefficient(pde.gravity,xq);
                if size(gravity,2) > 1
                    gravity = -gravity(:,end);
                end
                fq = -rho(:).*gravity(:).*gradS(:);
            end
            for a = 1:6
                F = F+accumarray(double(elem2dof(:,a)),...
                    w(q)*area.*fq(:)*phi(a),[Ndof,1]);
            end
        end
    end

    function [res,Ares,bedCoefficientRes,etaMinRes,etaMaxRes] = ...
            evaluateresidual(uk)
        [Kres,etaMinRes,etaMaxRes] = assembleviscous(uk);
        [Kbres,bedCoefficientRes] = assemblebed(uk);
        Fres = assembleforce;
        Ares = Kres+Kbres;
        r = Pperiodic'*(Ares*uk-Fres);
        FredRes = Pperiodic'*Fres;
        res = norm(r)/max(1,norm(FredRes));
    end

    function xq = evalpoint(lambda,elemIdx)
        if nargin < 2, elemIdx = (1:NT)'; end
        xq = lambda(1)*node(elem(elemIdx,1),:) + ...
             lambda(2)*node(elem(elemIdx,2),:) + ...
             lambda(3)*node(elem(elemIdx,3),:);
    end

    function edgeLength = localedgelength(elemIdx,e)
        localEdge = [2 3; 3 1; 1 2];
        v1 = node(elem(elemIdx,localEdge(e,1)),:);
        v2 = node(elem(elemIdx,localEdge(e,2)),:);
        edgeLength = sqrt(sum((v1-v2).^2,2));
    end

    function tangentX = localedgetangentx(elemIdx,e)
        localEdge = [2 3; 3 1; 1 2];
        v1 = node(elem(elemIdx,localEdge(e,1)),:);
        v2 = node(elem(elemIdx,localEdge(e,2)),:);
        tangent = v2-v1;
        tangentX = tangent(:,1)./sqrt(sum(tangent.^2,2));
    end

    function [wdof,pdof] = diagnosticwp(uk)
        gradU = recoverdofgradient(uk);
        ux = gradU(:,1);
        uz = gradU(:,2);
        q = dofNode(:,2)+option.slope*dofNode(:,1);
        wdof = diagnosticw(uk,ux,q);

        epsII = ux.^2+0.25*uz.^2;
        Acoef = coefficient(pde.A,dofNode);
        ncoef = coefficient(pde.n,dofNode);
        eta = 0.5.*Acoef.^(-1./ncoef).*...
            (epsII+option.eps_reg^2).^((1-ncoef)./(2*ncoef));
        rho = coefficient(pde.rho,dofNode);
        gravity = coefficient(pde.gravity,dofNode);
        if size(gravity,2) > 1
            gravity = -gravity(:,end);
        end
        depth = max(0,option.H-q);
        tauzz = -2*eta.*ux;
        pdof = rho(:).*gravity(:).*depth(:)+tauzz(:);
    end

    function gradDof = recoverdofgradient(uk)
        localLambda = [1 0 0; 0 1 0; 0 0 1; ...
                       0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
        gradSum = zeros(Ndof,2);
        weightSum = zeros(Ndof,1);
        for a = 1:6
            Dphi = p2gradient(localLambda(a,:),Dlambda);
            gradLocal = evalgrad(uk,elem2dof,Dphi);
            dof = double(elem2dof(:,a));
            gradSum(:,1) = gradSum(:,1)+accumarray(dof,area.*gradLocal(:,1),[Ndof,1]);
            gradSum(:,2) = gradSum(:,2)+accumarray(dof,area.*gradLocal(:,2),[Ndof,1]);
            weightSum = weightSum+accumarray(dof,area,[Ndof,1]);
        end
        gradDof = gradSum./max(weightSum,eps);
    end

    function wdof = diagnosticw(uk,ux,q)
        wdof = zeros(Ndof,1);
        tol = option.periodic_tol;
        xKey = round(dofNode(:,1)/tol);
        [~,~,xGroup] = unique(xKey);
        for group = 1:max(xGroup)
            idx = find(xGroup == group);
            [qUnique,~,qGroup] = unique(round(q(idx)/tol));
            qValue = qUnique*tol;
            uValue = accumarray(qGroup,uk(idx),[],@mean);
            uxValue = accumarray(qGroup,ux(idx),[],@mean);
            [qValue,order] = sort(qValue);
            uValue = uValue(order);
            uxValue = uxValue(order);
            wValue = zeros(size(qValue));
            wValue(1) = -option.slope*uValue(1);
            for iq = 2:length(qValue)
                dz = qValue(iq)-qValue(iq-1);
                wValue(iq) = wValue(iq-1)-0.5*dz*(uxValue(iq-1)+uxValue(iq));
            end
            inverseOrder = zeros(size(order));
            inverseOrder(order) = 1:length(order);
            wdof(idx) = wValue(inverseOrder(qGroup));
        end
    end
end

function [node,elem] = rectmesh2(box,h)
if isscalar(h), h = [h,h]; end
hx = h(1);
hz = h(2);
x = linspace(box(1),box(2),max(1,ceil((box(2)-box(1))/hx))+1);
z = linspace(box(3),box(4),max(1,ceil((box(4)-box(3))/hz))+1);
[xx,zz] = meshgrid(x,z);
node = [xx(:),zz(:)];
nz = length(z);
nx = length(x);
elem = zeros(2*(nx-1)*(nz-1),3);
idx = 1;
nodeidx = @(i,j) (i-1)*nz+j;
for i = 1:nx-1
    for j = 1:nz-1
        n1 = nodeidx(i,j);
        n2 = nodeidx(i+1,j);
        n3 = nodeidx(i+1,j+1);
        n4 = nodeidx(i,j+1);
        elem(idx:idx+1,:) = [n2 n3 n1; n4 n1 n3];
        idx = idx+2;
    end
end
end

function Dphi = p2gradient(l,Dlambda)
NT = size(Dlambda,1);
Dphi = zeros(NT,2,6);
Dphi(:,:,1) = (4*l(1)-1).*Dlambda(:,:,1);
Dphi(:,:,2) = (4*l(2)-1).*Dlambda(:,:,2);
Dphi(:,:,3) = (4*l(3)-1).*Dlambda(:,:,3);
Dphi(:,:,4) = 4*(l(2)*Dlambda(:,:,3)+l(3)*Dlambda(:,:,2));
Dphi(:,:,5) = 4*(l(3)*Dlambda(:,:,1)+l(1)*Dlambda(:,:,3));
Dphi(:,:,6) = 4*(l(1)*Dlambda(:,:,2)+l(2)*Dlambda(:,:,1));
end

function phi = p2basis(l)
phi = [l(1)*(2*l(1)-1), l(2)*(2*l(2)-1), ...
       l(3)*(2*l(3)-1), 4*l(2)*l(3), ...
       4*l(3)*l(1), 4*l(1)*l(2)];
end

function gradValue = evalgrad(value,elem2dof,Dphi)
NT = size(Dphi,1);
localValue = reshape(value(elem2dof),NT,1,6);
gradValue = zeros(NT,2);
for a = 1:6
    gradValue = gradValue+Dphi(:,:,a).*localValue(:,:,a);
end
end

function dofs = edgelocaldof(e)
switch e
    case 1
        dofs = [2 3 4];
    case 2
        dofs = [3 1 5];
    case 3
        dofs = [1 2 6];
end
end

function l = edgelambda(e,s)
switch e
    case 1
        l = [0,1-s,s];
    case 2
        l = [s,0,1-s];
    case 3
        l = [1-s,s,0];
end
end

function [s,w] = gaussedge(order)
if order <= 2
    s = [0.5-sqrt(1/12); 0.5+sqrt(1/12)];
    w = [0.5; 0.5];
else
    s = [0.5-sqrt(3/20); 0.5; 0.5+sqrt(3/20)];
    w = [5/18; 4/9; 5/18];
end
end

function [rep,info] = periodicrepresentatives(dofNode,periodicX,slope,tol)
Ndof = size(dofNode,1);
xl = periodicX(1);
xr = periodicX(2);
key = [dofNode(:,1), dofNode(:,2)+slope*dofNode(:,1)];
isRight = abs(key(:,1)-xr) <= tol*max(1,xr-xl);
key(isRight,1) = xl;
key = round(key/tol);
[~,ia,ic] = unique(key,'rows');
rep = ia(ic);
info.enabled = true;
info.direction = 1;
info.periodic_x = periodicX;
info.slope = slope;
info.masterDof = unique(rep);
info.numDof = Ndof;
info.numMasterDof = length(info.masterDof);
end

function value = coefficient(data,x)
if isa(data,'function_handle')
    value = data(x);
else
    value = data;
end
if isscalar(value)
    value = repmat(value,size(x,1),1);
elseif isrow(value) && size(x,1)>1
    value = repmat(value,size(x,1),1);
end
end

function value = getoption(option,name,defaultValue)
if isfield(option,name)
    value = option.(name);
else
    value = defaultValue;
end
end

function option = setoption(option,name,value)
if ~isfield(option,name), option.(name) = value; end
end
