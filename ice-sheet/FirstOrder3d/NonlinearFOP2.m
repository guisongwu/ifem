function [soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option)
%% NONLINEARFOP2 First-order Blatter-Pattyn ice-sheet model.
%
% Solve the three-dimensional first-order (FO/BP) ice-flow equations for
% the horizontal velocity U = (u,v) with quadratic tetrahedral elements.
% The weak form is assembled for -div(stress(U)) = f.  If pde.f is absent,
% the driving force is f = -rho*g*gradS, which corresponds to the strong
% form div(stress(U)) = rho*g*gradS.
%
% Boundary flags follow setboundary3:
%   1: Dirichlet horizontal velocity, pde.g_D defaults to zero.
%   2: prescribed horizontal traction, pde.g_N defaults to zero.
%   3: bed.  The default is Weertman sliding
%          beta*|U|^(m-1)*U
%      Set option.bed_condition = 'no-slip' to impose U = 0 instead.
%
% Required/optional PDE data:
%   pde.A, pde.n, pde.rho, pde.gravity, pde.gradS, pde.beta, pde.m,
%   pde.f, pde.g_D, pde.g_N.
%
% Important options:
%   option.tol, option.maxIt, option.damping, option.eps_reg,
%   option.quadorder, option.printlevel, option.bed_condition,
%   option.periodic, option.periodicBox, option.periodicSlope,
%   option.assemble_tangent.

if nargin < 5, option = struct; end
option = setoption(option,'tol',1e-8);
option = setoption(option,'maxIt',50);
option = setoption(option,'damping',0.7);
option = setoption(option,'eps_reg',1e-8);
option = setoption(option,'quadorder',4);
option = setoption(option,'facequadorder',4);
option = setoption(option,'printlevel',1);
option = setoption(option,'bed_condition','sliding');
option = setoption(option,'periodic',[]);
option = setoption(option,'periodicBox',[]);
option = setoption(option,'periodicSlope',[0,0]);
option = setoption(option,'periodic_tol',1e-10);
option = setoption(option,'assemble_tangent',false);
option = setoption(option,'residual_tol',option.tol);
option = setoption(option,'residual_check_threshold',...
    max(1e-2,sqrt(option.residual_tol)));

if ~isfield(pde,'A'), pde.A = 1; end
if ~isfield(pde,'n'), pde.n = 3; end
if ~isfield(pde,'rho'), pde.rho = 1; end
if ~isfield(pde,'gravity'), pde.gravity = 1; end
if ~isfield(pde,'beta'), pde.beta = 1; end
if ~isfield(pde,'m'), pde.m = 1; end
if ~isfield(pde,'g_N'), pde.g_N = []; end
if ~isfield(pde,'g_D'), pde.g_D = []; end
if ~isfield(pde,'gradS'), pde.gradS = []; end

if ~(ischar(option.bed_condition) || ...
        (isstring(option.bed_condition) && isscalar(option.bed_condition)))
    error('iFEM:FOBedCondition',...
        'bed_condition must be a character vector or scalar string.');
end
bedCondition = lower(strtrim(char(option.bed_condition)));
if ~ismember(bedCondition,{'sliding','no-slip'})
    error('iFEM:FOBedCondition',...
        'Unknown bed_condition value: %s.',bedCondition);
end

elem = double(elem);
[elem2dof,edge] = dof3P2(elem);
[Dlambda,volume] = gradbasis3(node,elem);
N = size(node,1);
NT = size(elem,1);
NE = size(edge,1);
Ndof = N+NE;
udofNode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
[periodicRep,periodicInfo] = periodicrepresentatives(udofNode,option);
[periodicMasterDof,~,scalarMaster] = unique(periodicRep);
NscalarMaster = length(periodicMasterDof);
state2master = [scalarMaster; NscalarMaster+scalarMaster];
Nmaster = 2*NscalarMaster;
Pperiodic = sparse((1:2*Ndof)',state2master,1,2*Ndof,Nmaster);

if isempty(bdFlag)
    bdFlag = zeros(NT,4,'uint8');
end

dirichletDof = finddirichletdof;
bedNoSlipDof = [];
if strcmp(bedCondition,'no-slip')
    bedNoSlipDof = findbeddof;
    dirichletDof = union(dirichletDof,bedNoSlipDof);
end
dirichletDof = dirichletDof(:);
fixedDof = [dirichletDof; Ndof+dirichletDof];
fixedMaster = unique(state2master(fixedDof));
freeMaster = setdiff((1:Nmaster)',fixedMaster);
freeDof = find(ismember(state2master,freeMaster));

state = zeros(2*Ndof,1);
state(fixedDof) = dirichletvalue(dirichletDof);
masterState = zeros(Nmaster,1);
if ~isempty(fixedDof)
    masterState(fixedMaster) = accumarray(state2master(fixedDof),...
        state(fixedDof),[Nmaster,1],@mean,0);
end
state = Pperiodic*masterState;
if isfield(option,'u0') && numel(option.u0) == 2*Ndof
    masterState = accumarray(state2master,option.u0(:),[Nmaster,1],@mean,0);
    state = Pperiodic*masterState;
    state(fixedDof) = dirichletvalue(dirichletDof);
    if ~isempty(fixedDof)
        masterState(fixedMaster) = accumarray(state2master(fixedDof),...
            state(fixedDof),[Nmaster,1],@mean,0);
    end
    state = Pperiodic*masterState;
end

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
    F = addtraction(F);
    A = A+Ab;
    Ared = Pperiodic'*A*Pperiodic;
    Fred = Pperiodic'*F;

    rhs = Fred-Ared(:,fixedMaster)*masterState(fixedMaster);
    newMasterState = masterState;
    newMasterState(freeMaster) = Ared(freeMaster,freeMaster)\rhs(freeMaster);
    newMasterState(fixedMaster) = masterState(fixedMaster);

    alpha = option.damping;
    updatedMasterState = masterState;
    updatedMasterState(freeMaster) = (1-alpha)*masterState(freeMaster)+...
        alpha*newMasterState(freeMaster);
    updatedMasterState(fixedMaster) = masterState(fixedMaster);
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
            fprintf(['FO P2 %2d: relchange %.3e, residual %.3e, ',...
                'eta [%.3e, %.3e]\n'],...
                k,residual(k),equationResidual(k),etaMin,etaMax);
        else
            fprintf('FO P2 %2d: relchange %.3e, eta [%.3e, %.3e]\n',...
                k,residual(k),etaMin,etaMax);
        end
    end
    if residualChecked(k) && residual(k) < option.tol && ...
            equationResidual(k) < option.residual_tol
        break
    end
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
info.bedCondition = bedCondition;
info.dirichletDof = dirichletDof;
info.periodic = periodicInfo;

soln.u = state(1:Ndof);
soln.v = state(Ndof+1:end);
soln.U = state;

eqn.A = A;
eqn.edge = edge;
eqn.elem2dof = elem2dof;
eqn.dofNode = udofNode;
eqn.freeDof = freeDof;
eqn.fixedDof = fixedDof;
eqn.periodicProjection = Pperiodic;
eqn.freeMaster = freeMaster;
eqn.fixedMaster = fixedMaster;
if option.assemble_tangent
    Kt = assembleviscoustangent(state);
    Kbt = assemblebedtangent(state);
    eqn.tangent = Pperiodic'*(Kt+Kbt)*Pperiodic;
    eqn.applyBetaDerivative = @assemblebetadirection;
end

    function [A,etaMin,etaMax] = assembleviscous(uk)
        [lambda,w] = quadpts3(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*400*NT,1);
        cols = zeros(nq*400*NT,1);
        vals = zeros(nq*400*NT,1);
        cursor = 0;
        etaMin = inf;
        etaMax = 0;
        uh = uk(1:Ndof);
        vh = uk(Ndof+1:end);
        localDof = [elem2dof,Ndof+elem2dof];

        for q = 1:nq
            Dphi = p2gradient3(lambda(q,:),Dlambda);
            gradU = evalgrad(uh,elem2dof,Dphi);
            gradV = evalgrad(vh,elem2dof,Dphi);
            ux = gradU(:,1); uy = gradU(:,2); uz = gradU(:,3);
            vx = gradV(:,1); vy = gradV(:,2); vz = gradV(:,3);
            epsII = ux.^2+vy.^2+ux.*vy+0.25*(uy+vx).^2+...
                0.25*uz.^2+0.25*vz.^2;
            xq = evalpoint(lambda(q,:));
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            eta = 0.5.*Aq.^(-1./nqfield).*...
                (epsII+option.eps_reg^2).^((1-nqfield)./(2*nqfield));
            etaMin = min(etaMin,min(eta));
            etaMax = max(etaMax,max(eta));

            for a = 1:20
                [aComp,aBasis] = splitlocal(a);
                da = Dphi(:,:,aBasis);
                for b = 1:20
                    [bComp,bBasis] = splitlocal(b);
                    db = Dphi(:,:,bBasis);
                    kab = w(q)*volume.*eta.*...
                        fostrainproduct(aComp,da,bComp,db);
                    idx = cursor+(1:NT);
                    rows(idx) = localDof(:,a);
                    cols(idx) = localDof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        A = sparse(rows,cols,vals,2*Ndof,2*Ndof);
    end

    function At = assembleviscoustangent(uk)
        [lambda,w] = quadpts3(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*400*NT,1);
        cols = zeros(nq*400*NT,1);
        vals = zeros(nq*400*NT,1);
        cursor = 0;
        uh = uk(1:Ndof);
        vh = uk(Ndof+1:end);
        localDof = [elem2dof,Ndof+elem2dof];

        for q = 1:nq
            Dphi = p2gradient3(lambda(q,:),Dlambda);
            gradU = evalgrad(uh,elem2dof,Dphi);
            gradV = evalgrad(vh,elem2dof,Dphi);
            ux = gradU(:,1); uy = gradU(:,2); uz = gradU(:,3);
            vx = gradV(:,1); vy = gradV(:,2); vz = gradV(:,3);
            epsII = ux.^2+vy.^2+ux.*vy+0.25*(uy+vx).^2+...
                0.25*uz.^2+0.25*vz.^2;
            xq = evalpoint(lambda(q,:));
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            strain = epsII+option.eps_reg^2;
            exponent = (1-nqfield)./(2*nqfield);
            eta = 0.5.*Aq.^(-1./nqfield).*strain.^exponent;

            for a = 1:20
                [aComp,aBasis] = splitlocal(a);
                da = Dphi(:,:,aBasis);
                stateTest = fostrainproduct(aComp,da,1,gradU)+...
                    fostrainproduct(aComp,da,2,gradV);
                for b = 1:20
                    [bComp,bBasis] = splitlocal(b);
                    db = Dphi(:,:,bBasis);
                    strainProduct = fostrainproduct(aComp,da,bComp,db);
                    stateDirection = epsdirection(bComp,db,ux,uy,uz,...
                        vx,vy,vz);
                    kab = w(q)*volume.*eta.*...
                        (strainProduct+exponent./strain.*...
                         stateTest.*stateDirection);
                    idx = cursor+(1:NT);
                    rows(idx) = localDof(:,a);
                    cols(idx) = localDof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        At = sparse(rows,cols,vals,2*Ndof,2*Ndof);
    end

    function [Ab,coefAtFace] = assemblebed(uk)
        Ab = sparse(2*Ndof,2*Ndof);
        coefAtFace = [];
        if strcmp(bedCondition,'no-slip') || ~any(bdFlag(:)==3)
            return;
        end
        [lambdaFace,wFace] = quadpts3face(option.facequadorder);
        nQuadPerFace = size(lambdaFace,1)/4;
        uh = uk(1:Ndof);
        vh = uk(Ndof+1:end);
        ii = [];
        jj = [];
        ss = [];
        coefAtFace = zeros(nnz(bdFlag(:)==3),1);
        coefCounter = 0;
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==3);
            if isempty(elemIdx), continue; end
            faceArea = localfacearea(elemIdx,f);
            for q = (f-1)*nQuadPerFace+(1:nQuadPerFace)
                phi = p2basis(lambdaFace(q,:));
                faceDof = elem2dof(elemIdx,:);
                xq = evalpoint(lambdaFace(q,:),elemIdx);
                beta = coefficient(pde.beta,xq);
                uq = uh(faceDof)*phi';
                vq = vh(faceDof)*phi';
                speed2 = uq.^2+vq.^2+option.eps_reg^2;
                gamma = beta.*speed2.^((pde.m-1)/2);
                if q == (f-1)*nQuadPerFace+ceil(nQuadPerFace/2)
                    coefAtFace(coefCounter+(1:length(elemIdx))) = gamma;
                    coefCounter = coefCounter+length(elemIdx);
                end
                for a = 1:10
                    ia = faceDof(:,a);
                    for b = 1:10
                        ib = faceDof(:,b);
                        s = wFace(q)*faceArea.*gamma*phi(a)*phi(b);
                        ii = [ii;ia;Ndof+ia]; %#ok<AGROW>
                        jj = [jj;ib;Ndof+ib]; %#ok<AGROW>
                        ss = [ss;s;s]; %#ok<AGROW>
                    end
                end
            end
        end
        Ab = sparse(double(ii),double(jj),ss,2*Ndof,2*Ndof);
    end

    function Abt = assemblebedtangent(uk)
        Abt = sparse(2*Ndof,2*Ndof);
        if strcmp(bedCondition,'no-slip') || ~any(bdFlag(:)==3)
            return;
        end
        [lambdaFace,wFace] = quadpts3face(option.facequadorder);
        nQuadPerFace = size(lambdaFace,1)/4;
        uh = uk(1:Ndof);
        vh = uk(Ndof+1:end);
        ii = [];
        jj = [];
        ss = [];
        exponent = (pde.m-1)/2;
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==3);
            if isempty(elemIdx), continue; end
            faceArea = localfacearea(elemIdx,f);
            for q = (f-1)*nQuadPerFace+(1:nQuadPerFace)
                phi = p2basis(lambdaFace(q,:));
                faceDof = elem2dof(elemIdx,:);
                xq = evalpoint(lambdaFace(q,:),elemIdx);
                beta = coefficient(pde.beta,xq);
                uq = uh(faceDof)*phi';
                vq = vh(faceDof)*phi';
                speed2 = uq.^2+vq.^2+option.eps_reg^2;
                base = speed2.^exponent;
                cross = (pde.m-1)*speed2.^(exponent-1);
                block11 = beta.*(base+cross.*uq.^2);
                block12 = beta.*cross.*uq.*vq;
                block22 = beta.*(base+cross.*vq.^2);
                for a = 1:10
                    ia = faceDof(:,a);
                    for b = 1:10
                        ib = faceDof(:,b);
                        mass = wFace(q)*faceArea*phi(a)*phi(b);
                        ii = [ii;ia;ia;Ndof+ia;Ndof+ia]; %#ok<AGROW>
                        jj = [jj;ib;Ndof+ib;ib;Ndof+ib]; %#ok<AGROW>
                        ss = [ss;mass.*block11;mass.*block12;...
                            mass.*block12;mass.*block22]; %#ok<AGROW>
                    end
                end
            end
        end
        Abt = sparse(double(ii),double(jj),ss,2*Ndof,2*Ndof);
    end

    function load = assemblebetadirection(betaDirection)
        fullLoad = zeros(2*Ndof,1);
        if strcmp(bedCondition,'no-slip') || ~any(bdFlag(:)==3)
            load = Pperiodic'*fullLoad;
            return;
        end
        [lambdaFace,wFace] = quadpts3face(option.facequadorder);
        nQuadPerFace = size(lambdaFace,1)/4;
        uh = state(1:Ndof);
        vh = state(Ndof+1:end);
        exponent = (pde.m-1)/2;
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==3);
            if isempty(elemIdx), continue; end
            faceArea = localfacearea(elemIdx,f);
            faceDof = elem2dof(elemIdx,:);
            for q = (f-1)*nQuadPerFace+(1:nQuadPerFace)
                phi = p2basis(lambdaFace(q,:));
                xq = evalpoint(lambdaFace(q,:),elemIdx);
                deltaBeta = coefficient(betaDirection,xq);
                uq = uh(faceDof)*phi';
                vq = vh(faceDof)*phi';
                speed = (uq.^2+vq.^2+option.eps_reg^2).^exponent;
                directionU = deltaBeta.*speed.*uq;
                directionV = deltaBeta.*speed.*vq;
                for a = 1:10
                    contributionU = wFace(q)*faceArea.*...
                        directionU*phi(a);
                    contributionV = wFace(q)*faceArea.*...
                        directionV*phi(a);
                    fullLoad = fullLoad+accumarray(...
                        double(faceDof(:,a)),contributionU,...
                        [2*Ndof,1]);
                    fullLoad = fullLoad+accumarray(...
                        double(Ndof+faceDof(:,a)),contributionV,...
                        [2*Ndof,1]);
                end
            end
        end
        load = Pperiodic'*fullLoad;
    end

    function F = assembleforce
        F = zeros(2*Ndof,1);
        [lambda,w] = quadpts3(option.quadorder);
        for q = 1:size(lambda,1)
            phi = p2basis(lambda(q,:));
            xq = evalpoint(lambda(q,:));
            if isfield(pde,'f') && ~isempty(pde.f)
                fq = coefficient(pde.f,xq);
                if size(fq,2) ~= 2
                    error('pde.f must return an N-by-2 array.');
                end
            else
                gradS = surfacedrivinggradient(xq);
                rho = coefficient(pde.rho,xq);
                gravity = coefficient(pde.gravity,xq);
                if size(gravity,2) > 1
                    gravity = gravity(:,end);
                end
                fq = -[rho(:).*gravity(:).*gradS(:,1),...
                       rho(:).*gravity(:).*gradS(:,2)];
            end
            for a = 1:10
                F = F+accumarray(double(elem2dof(:,a)),...
                    w(q)*volume.*fq(:,1)*phi(a),[2*Ndof,1]);
                F = F+accumarray(double(Ndof+elem2dof(:,a)),...
                    w(q)*volume.*fq(:,2)*phi(a),[2*Ndof,1]);
            end
        end
    end

    function F = addtraction(F)
        if isempty(pde.g_N) || ~any(bdFlag(:)==2), return; end
        [lambdaFace,wFace] = quadpts3face(option.facequadorder);
        nQuadPerFace = size(lambdaFace,1)/4;
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==2);
            if isempty(elemIdx), continue; end
            faceArea = localfacearea(elemIdx,f);
            faceDof = elem2dof(elemIdx,:);
            for q = (f-1)*nQuadPerFace+(1:nQuadPerFace)
                phi = p2basis(lambdaFace(q,:));
                xq = evalpoint(lambdaFace(q,:),elemIdx);
                tq = coefficient(pde.g_N,xq);
                if size(tq,2) ~= 2
                    error('pde.g_N must return an N-by-2 array.');
                end
                for a = 1:10
                    F = F+accumarray(double(faceDof(:,a)),...
                        wFace(q)*faceArea.*tq(:,1)*phi(a),[2*Ndof,1]);
                    F = F+accumarray(double(Ndof+faceDof(:,a)),...
                        wFace(q)*faceArea.*tq(:,2)*phi(a),[2*Ndof,1]);
                end
            end
        end
    end

    function [res,Ares,bedCoefficientRes,etaMinRes,etaMaxRes] = ...
            evaluateresidual(uk)
        [Kres,etaMinRes,etaMaxRes] = assembleviscous(uk);
        [Kbres,bedCoefficientRes] = assemblebed(uk);
        Fres = assembleforce;
        Fres = addtraction(Fres);
        Ares = Kres+Kbres;
        r = Pperiodic'*(Ares*uk-Fres);
        FredRes = Pperiodic'*Fres;
        res = norm(r(freeMaster))/max(1,norm(FredRes(freeMaster)));
    end

    function dof = finddirichletdof
        dof = [];
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==1);
            if isempty(elemIdx), continue; end
            localDof = facelocaldof(f);
            dof = [dof; elem2dof(elemIdx,localDof)']; %#ok<AGROW>
        end
        dof = unique(double(dof(:)));
    end

    function dof = findbeddof
        dof = [];
        for f = 1:4
            elemIdx = find(bdFlag(:,f)==3);
            if isempty(elemIdx), continue; end
            localDof = facelocaldof(f);
            dof = [dof; elem2dof(elemIdx,localDof)']; %#ok<AGROW>
        end
        dof = unique(double(dof(:)));
    end

    function value = dirichletvalue(dof)
        value = zeros(2*length(dof),1);
        if isempty(dof) || isempty(pde.g_D), return; end
        g = coefficient(pde.g_D,udofNode(dof,:));
        if size(g,2) ~= 2
            error('pde.g_D must return an N-by-2 array.');
        end
        if ~isempty(bedNoSlipDof)
            g(ismember(dof,bedNoSlipDof),:) = 0;
        end
        value = [g(:,1);g(:,2)];
    end

    function xq = evalpoint(lambda,elemIdx)
        if nargin < 2, elemIdx = (1:NT)'; end
        xq = lambda(1)*node(elem(elemIdx,1),:) + ...
             lambda(2)*node(elem(elemIdx,2),:) + ...
             lambda(3)*node(elem(elemIdx,3),:) + ...
             lambda(4)*node(elem(elemIdx,4),:);
    end

    function gradS = surfacedrivinggradient(xq)
        if isfield(pde,'gradS') && ~isempty(pde.gradS)
            gradS = coefficient(pde.gradS,xq);
        elseif isfield(pde,'s') && ~isempty(pde.s)
            h = sqrt(eps)*max(1,max(abs(xq(:))));
            xp = xq; xm = xq;
            xp(:,1) = xp(:,1)+h; xm(:,1) = xm(:,1)-h;
            sx = (coefficient(pde.s,xp)-coefficient(pde.s,xm))/(2*h);
            xp = xq; xm = xq;
            xp(:,2) = xp(:,2)+h; xm(:,2) = xm(:,2)-h;
            sy = (coefficient(pde.s,xp)-coefficient(pde.s,xm))/(2*h);
            gradS = [sx(:),sy(:)];
        else
            gradS = zeros(size(xq,1),2);
        end
        if size(gradS,2) ~= 2
            error('pde.gradS must return an N-by-2 array.');
        end
    end

    function area = localfacearea(elemIdx,f)
        face = localface(f);
        p1 = node(elem(elemIdx,face(1)),:);
        p2 = node(elem(elemIdx,face(2)),:);
        p3 = node(elem(elemIdx,face(3)),:);
        area = sqrt(sum(mycross(p2-p1,p3-p1,2).^2,2))/2;
    end
end

function Dphi = p2gradient3(l,Dlambda)
NT = size(Dlambda,1);
Dphi = zeros(NT,3,10);
Dphi(:,:,1) = (4*l(1)-1).*Dlambda(:,:,1);
Dphi(:,:,2) = (4*l(2)-1).*Dlambda(:,:,2);
Dphi(:,:,3) = (4*l(3)-1).*Dlambda(:,:,3);
Dphi(:,:,4) = (4*l(4)-1).*Dlambda(:,:,4);
Dphi(:,:,5) = 4*(l(1)*Dlambda(:,:,2)+l(2)*Dlambda(:,:,1));
Dphi(:,:,6) = 4*(l(1)*Dlambda(:,:,3)+l(3)*Dlambda(:,:,1));
Dphi(:,:,7) = 4*(l(1)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,1));
Dphi(:,:,8) = 4*(l(2)*Dlambda(:,:,3)+l(3)*Dlambda(:,:,2));
Dphi(:,:,9) = 4*(l(2)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,2));
Dphi(:,:,10) = 4*(l(3)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,3));
end

function phi = p2basis(l)
phi = [l(1)*(2*l(1)-1), l(2)*(2*l(2)-1), ...
       l(3)*(2*l(3)-1), l(4)*(2*l(4)-1), ...
       4*l(1)*l(2), 4*l(1)*l(3), 4*l(1)*l(4), ...
       4*l(2)*l(3), 4*l(2)*l(4), 4*l(3)*l(4)];
end

function gradValue = evalgrad(value,elem2dof,Dphi)
NT = size(Dphi,1);
localValue = reshape(value(elem2dof),NT,1,10);
gradValue = zeros(NT,3);
for a = 1:10
    gradValue = gradValue+Dphi(:,:,a).*localValue(:,:,a);
end
end

function product = fostrainproduct(testComp,da,trialComp,db)
if testComp == 1 && trialComp == 1
    product = 4*da(:,1).*db(:,1)+da(:,2).*db(:,2)+...
        da(:,3).*db(:,3);
elseif testComp == 1 && trialComp == 2
    product = 2*da(:,1).*db(:,2)+da(:,2).*db(:,1);
elseif testComp == 2 && trialComp == 1
    product = da(:,1).*db(:,2)+2*da(:,2).*db(:,1);
else
    product = da(:,1).*db(:,1)+4*da(:,2).*db(:,2)+...
        da(:,3).*db(:,3);
end
end

function value = epsdirection(component,db,ux,uy,uz,vx,vy,vz)
if component == 1
    value = (2*ux+vy).*db(:,1)+0.5*(uy+vx).*db(:,2)+...
        0.5*uz.*db(:,3);
else
    value = 0.5*(uy+vx).*db(:,1)+(ux+2*vy).*db(:,2)+...
        0.5*vz.*db(:,3);
end
end

function [component,basis] = splitlocal(a)
component = 1+(a>10);
basis = a-10*(component-1);
end

function dofs = facelocaldof(f)
switch f
    case 1
        dofs = [2 3 4 8 9 10];
    case 2
        dofs = [1 3 4 6 7 10];
    case 3
        dofs = [1 2 4 5 7 9];
    case 4
        dofs = [1 2 3 5 6 8];
end
end

function face = localface(f)
switch f
    case 1
        face = [2 4 3];
    case 2
        face = [1 3 4];
    case 3
        face = [1 4 2];
    case 4
        face = [1 2 3];
end
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

function option = setoption(option,name,value)
if ~isfield(option,name), option.(name) = value; end
end

function [rep,info] = periodicrepresentatives(dofNode,option)
Ndof = size(dofNode,1);
rep = (1:Ndof)';
info.enabled = false;
info.directions = [];
info.masterDof = rep;
if ~isfield(option,'periodic') || isempty(option.periodic) || ...
        (islogical(option.periodic) && ~any(option.periodic))
    return;
end
if isempty(option.periodicBox)
    error('iFEM:FOPeriodicBox',...
        'option.periodicBox is required when option.periodic is set.');
end
periodic = option.periodic;
if islogical(periodic)
    periodic = find(periodic);
end
box = option.periodicBox;
if isvector(box)
    box = reshape(box,2,[])';
end
if size(box,1) == 2 && size(box,2) == 3
    box = box';
end
if size(box,1) ~= 3 || size(box,2) ~= 2
    error('iFEM:FOPeriodicBox',...
        'option.periodicBox must be a 3-by-2 array or [x0 x1 y0 y1 z0 z1].');
end
tol = option.periodic_tol;
key = dofNode;
periodicSlope = option.periodicSlope;
if isscalar(periodicSlope)
    periodicSlope = [periodicSlope,0];
end
for d = periodic(:)'
    lower = box(d,1);
    upper = box(d,2);
    isUpper = abs(key(:,d)-upper) <= tol*max(1,abs(upper-lower));
    if d <= 2
        key(isUpper,3) = key(isUpper,3)+...
            periodicSlope(d)*(upper-lower);
    end
    key(isUpper,d) = lower;
end
key = round(key/tol);
[~,ia,ic] = unique(key,'rows');
rep = ia(ic);
info.enabled = true;
info.directions = periodic(:)';
info.masterDof = unique(rep);
info.numDof = Ndof;
info.numMasterDof = length(info.masterDof);
end
