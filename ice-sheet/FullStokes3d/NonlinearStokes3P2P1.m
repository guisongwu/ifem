function [soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option)
%% NONLINEARSTOKES3P2P1 Nonlinear full-Stokes model on a periodic cuboid.
%
% Solve
%
%   div(u) = 0,
%   -div(eta(u)*(grad(u)+grad(u)')) + grad(p) = f,
%
% with P2 velocity and P1 pressure on tetrahedra in a rectangular periodic
% slab.  The x=periodic_x(1)/periodic_x(2) and
% y=periodic_y(1)/periodic_y(2) faces are identified by default.  For a
% tilted slab
%
%   z = z0 - sx*x - sy*y,
%
% set option.periodic_slope = [sx,sy] (or the legacy fields
% option.periodic_slope_x and option.periodic_slope_y).  Periodic matching
% is then performed in the coordinates (x,y,z+sx*x+sy*y).
%
% Boundary flags are
%
%   1: no-slip Dirichlet u = 0,
%   2: prescribed traction sigma*n = g_N,
%   3: impermeable Weertman bed, u*n = 0 and
%      beta*|u_t|^(m-1)*u_t + tangential traction = 0.
%
% Required data are pde.A, pde.n, pde.beta, and pde.m.  Optional data are
% pde.f, pde.rho, pde.gravity, pde.g_N, pde.beta_scale, and
% pde.pressure_dof_scale.  Set option.assemble_tangent true when the
% consistent nonlinear Jacobian is needed by an adjoint solve.
% Pressure is determined up to a constant.  Use
% option.pressure_constraint = 'pin' to fix one pressure degree of freedom
% during the linear solve; the returned pressure is shifted to zero mean.
% Do not add pressure constraints as Lagrange multipliers, since they enter
% the discrete divergence equation.

if nargin < 5, option = struct; end
option = setoption(option,'tol',1e-8);
option = setoption(option,'maxIt',50);
option = setoption(option,'damping',0.7);
option = setoption(option,'eps_reg',1e-8);
option = setoption(option,'quadorder',4);
option = setoption(option,'facequadorder',4);
option = setoption(option,'printlevel',1);
option = setoption(option,'assemble_tangent',false);
option = setoption(option,'pressure_constraint','none');
option = setoption(option,'periodic',true);
option = setoption(option,'periodic_x',[min(node(:,1)),max(node(:,1))]);
option = setoption(option,'periodic_y',[min(node(:,2)),max(node(:,2))]);
option = setoption(option,'periodic_slope_x',0);
option = setoption(option,'periodic_slope_y',0);
if isfield(option,'periodic_slope')
    if numel(option.periodic_slope) ~= 2
        error('iFEM:NS3PeriodicSlope',...
            'option.periodic_slope must contain [slope_x,slope_y].');
    end
    option.periodic_slope_x = option.periodic_slope(1);
    option.periodic_slope_y = option.periodic_slope(2);
else
    option.periodic_slope = [option.periodic_slope_x,option.periodic_slope_y];
end
option = setoption(option,'residual_tol',option.tol);
option = setoption(option,'residual_check_threshold',...
    max(1e-2,sqrt(option.residual_tol)));

if ~isfield(pde,'A'), pde.A = 1; end
if ~isfield(pde,'n'), pde.n = 3; end
if ~isfield(pde,'beta'), pde.beta = 1; end
if ~isfield(pde,'m'), pde.m = 1; end
if ~isfield(pde,'g_N'), pde.g_N = []; end
betaScale = getfielddefault(pde,'beta_scale',1);
pressureDofScale = getfielddefault(pde,'pressure_dof_scale',1);

[elem2dof,edge] = dof3P2(elem);
[Dlambda,volume] = gradbasis3(node,elem);
N = size(node,1);
NT = size(elem,1);
NE = size(edge,1);
Nu = N+NE;
Np = N;
udofNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];

isBedFace = bdFlag==3;
isTopFace = bdFlag==2;
isDirichletFace = bdFlag==1;
bedFace = find(isBedFace(:));
topFace = find(isTopFace(:));

dirichletDof = boundarydofs(isDirichletFace);
bedDof = boundarydofs(isBedFace);
bedNormalAtDof = beddofnormal(bedDof);

pressureConstraint = lower(strtrim(char(option.pressure_constraint)));
if ~ismember(pressureConstraint,{'none','pin','mean-zero'})
    error('iFEM:NS3PressureConstraint',...
        'pressure_constraint must be pin, mean-zero, or none.');
end
if strcmp(pressureConstraint,'mean-zero')
    error('iFEM:NS3PressureConstraint',...
        ['mean-zero pressure constraints are not supported in ',...
         'NonlinearStokes3P2P1 because multiplier pressure ',...
         'constraints pollute the divergence equation. Use ',...
         'pressure_constraint=''pin'' to remove the pressure null space.']);
end
addPressureMeanConstraint = false;

[C,bedNormal] = buildconstraints;
nConstraint = size(C,1);
B = Bmatrix;

u = zeros(3*Nu,1);
p = zeros(Np,1);
constraintMultiplier = zeros(nConstraint,1);
if isfield(option,'u0') && numel(option.u0) == 3*Nu
    u = option.u0(:);
end

residual = zeros(option.maxIt,1);
viscosityRange = zeros(option.maxIt,2);
momentumResidual = NaN(option.maxIt,1);
divergenceResidual = NaN(option.maxIt,1);
constraintResidual = NaN(option.maxIt,1);
nonlinearResidual = NaN(option.maxIt,1);
residualChecked = false(option.maxIt,1);
checkResidualEveryStep = false;
t0 = cputime;

for k = 1:option.maxIt
    [K,etaMin,etaMax] = assembleviscous(u);
    [Kb,bedCoefficient] = assemblebed(u,false);
    F = assembleforce;
    F = addtraction(F);

    M = [K+Kb,B';B,sparse(Np,Np)];
    rhs = [F;zeros(Np,1)];
    saddle = [M,C';C,sparse(nConstraint,nConstraint)];
    rhsFull = [rhs;zeros(nConstraint,1)];
    if strcmp(pressureConstraint,'pin')
        [saddle,rhsFull] = pinpressure(saddle,rhsFull);
    end
    fullsol = saddle\rhsFull;

    unew = fullsol(1:3*Nu);
    pnew = fullsol(3*Nu+(1:Np));
    multiplierNew = fullsol(3*Nu+Np+(1:nConstraint));
    alpha = option.damping;
    updatedU = (1-alpha)*u+alpha*unew;
    updatedP = (1-alpha)*p+alpha*pnew;
    updatedMultiplier = (1-alpha)*constraintMultiplier+alpha*multiplierNew;

    residual(k) = norm(updatedU-u)/max(1,norm(updatedU));
    viscosityRange(k,:) = [etaMin,etaMax];
    u = updatedU;
    p = updatedP;
    constraintMultiplier = updatedMultiplier;

    if residual(k) <= option.residual_check_threshold
        checkResidualEveryStep = true;
    end
    if checkResidualEveryStep || k == option.maxIt
        [momentumResidual(k),divergenceResidual(k),...
            constraintResidual(k),nonlinearResidual(k),K,Kb,...
            etaMin,etaMax,bedCoefficient] = ...
            evaluateresidual(u,p,constraintMultiplier);
        viscosityRange(k,:) = [etaMin,etaMax];
        residualChecked(k) = true;
    end

    if option.printlevel > 0
        if residualChecked(k)
            fprintf(['nonlinear Stokes3 %2d: relchange %.3e, ',...
                'residual %.3e, eta [%.3e, %.3e]\n'],...
                k,residual(k),nonlinearResidual(k),etaMin,etaMax);
        else
            fprintf('nonlinear Stokes3 %2d: relchange %.3e, eta [%.3e, %.3e]\n',...
                k,residual(k),etaMin,etaMax);
        end
    end
    if residualChecked(k) && residual(k) < option.tol && ...
            nonlinearResidual(k) < option.residual_tol
        break
    end
end

info.itStep = k;
info.relchange = residual(1:k);
info.viscosityRange = viscosityRange(1:k,:);
info.momentumResidual = momentumResidual(1:k);
info.divergenceResidual = divergenceResidual(1:k);
info.constraintResidual = constraintResidual(1:k);
info.nonlinearResidual = nonlinearResidual(1:k);
info.residualChecked = residualChecked(1:k);
info.residualTolerance = option.residual_tol;
info.converged = residualChecked(k) && residual(k) < option.tol && ...
    nonlinearResidual(k) < option.residual_tol;
info.solveTime = cputime-t0;
info.bedCoefficient = bedCoefficient;
info.pressureConstraint = pressureConstraint;
info.periodic = option.periodic;
info.periodicSlope = option.periodic_slope;

soln.u = u;
if strcmp(pressureConstraint,'pin')
    p = p-mean(p);
end
soln.p = p;
soln.p_phgism = p/pressureDofScale;
soln.ux = u(1:Nu);
soln.uy = u(Nu+(1:Nu));
soln.uz = u(2*Nu+(1:Nu));

eqn.K = K+Kb;
eqn.B = B;
eqn.C = C;
eqn.edge = edge;
eqn.elem2dof = elem2dof;
eqn.bedDof = bedDof;
eqn.bedNormal = bedNormal;
eqn.constraintMultiplier = constraintMultiplier;
eqn.udofNode = udofNode;

if option.assemble_tangent
    Kt = assembleviscoustangent(u);
    Kbt = assemblebed(u,true);
    tangentM = [Kt+Kbt,B';B,sparse(Np,Np)];
    eqn.tangent = [tangentM,C';C,sparse(nConstraint,nConstraint)];
    if strcmp(pressureConstraint,'pin')
        eqn.tangent = pinpressure(eqn.tangent,[]);
    end
    eqn.applyBetaDerivative = @assemblebetadirection;
end

    function [A,b] = pinpressure(A,b)
        pressureDof = 3*Nu+1;
        A(pressureDof,:) = 0;
        A(:,pressureDof) = 0;
        A(pressureDof,pressureDof) = 1;
        if ~isempty(b)
            b(pressureDof) = 0;
        end
    end

    function [K,etaMin,etaMax] = assembleviscous(uk)
        [lambda,w] = quadpts3(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*900*NT,1);
        cols = zeros(nq*900*NT,1);
        vals = zeros(nq*900*NT,1);
        cursor = 0;
        etaMin = inf;
        etaMax = 0;
        ux = uk(1:Nu);
        uy = uk(Nu+(1:Nu));
        uz = uk(2*Nu+(1:Nu));
        for q = 1:nq
            Dphi = p2gradient3(lambda(q,:),Dlambda);
            [exx,eyy,ezz,exy,exz,eyz] = strainat(ux,uy,uz,Dphi);
            epsII = 0.5*(exx.^2+eyy.^2+ezz.^2+...
                2*(exy.^2+exz.^2+eyz.^2));
            xq = barycenter(lambda(q,:),node,elem);
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            eta = 0.5.*Aq.^(-1./nqfield).*...
                (epsII+option.eps_reg^2).^((1-nqfield)./(2*nqfield));
            etaMin = min(etaMin,min(eta));
            etaMax = max(etaMax,max(eta));
            localDof = [elem2dof,Nu+elem2dof,2*Nu+elem2dof];
            for a = 1:30
                [aComp,aBasis] = splitlocal(a);
                aStrain = basisstrain3(aComp,Dphi,aBasis);
                for b = 1:30
                    [bComp,bBasis] = splitlocal(b);
                    bStrain = basisstrain3(bComp,Dphi,bBasis);
                    strainDot = straininner(aStrain,bStrain);
                    kab = 2*w(q)*volume.*eta.*strainDot;
                    idx = cursor+(1:NT);
                    rows(idx) = localDof(:,a);
                    cols(idx) = localDof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        K = sparse(rows,cols,vals,3*Nu,3*Nu);
    end

    function Kt = assembleviscoustangent(uk)
        [lambda,w] = quadpts3(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*900*NT,1);
        cols = zeros(nq*900*NT,1);
        vals = zeros(nq*900*NT,1);
        cursor = 0;
        ux = uk(1:Nu);
        uy = uk(Nu+(1:Nu));
        uz = uk(2*Nu+(1:Nu));
        for q = 1:nq
            Dphi = p2gradient3(lambda(q,:),Dlambda);
            stateStrain = strainat(ux,uy,uz,Dphi);
            epsII = 0.5*(stateStrain{1}.^2+stateStrain{2}.^2+...
                stateStrain{3}.^2+2*(stateStrain{4}.^2+...
                stateStrain{5}.^2+stateStrain{6}.^2));
            xq = barycenter(lambda(q,:),node,elem);
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            strainRegularized = epsII+option.eps_reg^2;
            exponent = (1-nqfield)./(2*nqfield);
            eta = 0.5.*Aq.^(-1./nqfield).*strainRegularized.^exponent;
            localDof = [elem2dof,Nu+elem2dof,2*Nu+elem2dof];
            for a = 1:30
                [aComp,aBasis] = splitlocal(a);
                aStrain = basisstrain3(aComp,Dphi,aBasis);
                stateDotA = straininner(stateStrain,aStrain);
                for b = 1:30
                    [bComp,bBasis] = splitlocal(b);
                    bStrain = basisstrain3(bComp,Dphi,bBasis);
                    stateDotB = straininner(stateStrain,bStrain);
                    strainDot = straininner(aStrain,bStrain);
                    kab = 2*w(q)*volume.*eta.*...
                        (strainDot+exponent./strainRegularized.*...
                         stateDotA.*stateDotB);
                    idx = cursor+(1:NT);
                    rows(idx) = localDof(:,a);
                    cols(idx) = localDof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        Kt = sparse(rows,cols,vals,3*Nu,3*Nu);
    end

    function B = Bmatrix
        [lambda,w] = quadpts3(2);
        Bx = sparse(Np,Nu);
        By = sparse(Np,Nu);
        Bz = sparse(Np,Nu);
        for q = 1:size(lambda,1)
            Dphi = p2gradient3(lambda(q,:),Dlambda);
            for a = 1:10
                for j = 1:4
                    v = w(q)*volume.*lambda(q,j);
                    Bx = Bx+sparse(double(elem(:,j)),double(elem2dof(:,a)),...
                        -v.*Dphi(:,1,a),Np,Nu);
                    By = By+sparse(double(elem(:,j)),double(elem2dof(:,a)),...
                        -v.*Dphi(:,2,a),Np,Nu);
                    Bz = Bz+sparse(double(elem(:,j)),double(elem2dof(:,a)),...
                        -v.*Dphi(:,3,a),Np,Nu);
                end
            end
        end
        B = [Bx,By,Bz];
    end

    function [Kb,coefAtMid] = assemblebed(uk,tangentFlag)
        Kb = sparse(3*Nu,3*Nu);
        coefAtMid = [];
        if isempty(bedFace), return; end
        [lambdaAll,wAll] = quadpts3face(option.facequadorder);
        faceLocal = ceil((1:size(lambdaAll,1))'/(size(lambdaAll,1)/4));
        exponent = (pde.m-1)/2;
        ii = [];
        jj = [];
        ss = [];
        for id = 1:numel(bedFace)
            [t,f] = ind2sub(size(bdFlag),bedFace(id));
            area = localfacearea(t,f);
            normal = localfacenormal(t,f);
            projector = eye(3)-normal(:)*normal(:)';
            localDof = double(elem2dof(t,:));
            ux = uk(localDof);
            uy = uk(Nu+localDof);
            uz = uk(2*Nu+localDof);
            for q = find(faceLocal==f)'
                phi = p2value3(lambdaAll(q,:));
                xq = lambdaAll(q,1)*node(elem(t,1),:) + ...
                     lambdaAll(q,2)*node(elem(t,2),:) + ...
                     lambdaAll(q,3)*node(elem(t,3),:) + ...
                     lambdaAll(q,4)*node(elem(t,4),:);
                beta = betaScale*coefficient(pde.beta,xq);
                uq = [ux'*phi';uy'*phi';uz'*phi'];
                ut = projector*uq;
                ut2 = sum(ut.^2);
                speedRegularized = ut2+option.eps_reg^2;
                if tangentFlag
                    bedTensor = beta.*(speedRegularized.^exponent.*projector + ...
                        (pde.m-1).*speedRegularized.^(exponent-1).*(ut*ut'));
                else
                    bedTensor = beta.*speedRegularized.^exponent.*projector;
                end
                if q == find(faceLocal==f,1,'first')
                    coefAtMid = [coefAtMid;beta.*speedRegularized.^exponent]; %#ok<AGROW>
                end
                for a = 1:10
                    ia = localDof(a);
                    for b = 1:10
                        ib = localDof(b);
                        s = wAll(q)*area*phi(a)*phi(b);
                        rowDof = [ia;Nu+ia;2*Nu+ia];
                        colDof = [ib;Nu+ib;2*Nu+ib];
                        [rr,cc] = ndgrid(rowDof,colDof);
                        ii = [ii;rr(:)]; %#ok<AGROW>
                        jj = [jj;cc(:)]; %#ok<AGROW>
                        ss = [ss;s*bedTensor(:)]; %#ok<AGROW>
                    end
                end
            end
        end
        Kb = sparse(double(ii),double(jj),ss,3*Nu,3*Nu);
    end

    function load = assemblebetadirection(betaDirection)
        load = zeros(3*Nu+Np+nConstraint,1);
        if isempty(bedFace), return; end
        [lambdaAll,wAll] = quadpts3face(option.facequadorder);
        faceLocal = ceil((1:size(lambdaAll,1))'/(size(lambdaAll,1)/4));
        uxAll = u(1:Nu);
        uyAll = u(Nu+(1:Nu));
        exponent = (pde.m-1)/2;
        for id = 1:numel(bedFace)
            [t,f] = ind2sub(size(bdFlag),bedFace(id));
            area = localfacearea(t,f);
            normal = localfacenormal(t,f);
            projector = eye(3)-normal(:)*normal(:)';
            localDof = double(elem2dof(t,:));
            ux = uxAll(localDof);
            uy = uyAll(localDof);
            uzAll = u(2*Nu+(1:Nu));
            uz = uzAll(localDof);
            for q = find(faceLocal==f)'
                phi = p2value3(lambdaAll(q,:));
                xq = lambdaAll(q,1)*node(elem(t,1),:) + ...
                     lambdaAll(q,2)*node(elem(t,2),:) + ...
                     lambdaAll(q,3)*node(elem(t,3),:) + ...
                     lambdaAll(q,4)*node(elem(t,4),:);
                deltaBeta = betaScale*coefficient(betaDirection,xq);
                uq = [ux'*phi';uy'*phi';uz'*phi'];
                ut = projector*uq;
                factor = deltaBeta.*(sum(ut.^2)+option.eps_reg^2).^exponent;
                for a = 1:10
                    contribution = wAll(q)*area*factor*phi(a);
                    ia = localDof(a);
                    load(ia) = load(ia)+contribution*ut(1);
                    load(Nu+ia) = load(Nu+ia)+contribution*ut(2);
                    load(2*Nu+ia) = load(2*Nu+ia)+contribution*ut(3);
                end
            end
        end
    end

    function F = assembleforce
        F = zeros(3*Nu,1);
        [lambda,w] = quadpts3(option.quadorder);
        for q = 1:size(lambda,1)
            phi = p2value3(lambda(q,:));
            xq = barycenter(lambda(q,:),node,elem);
            if isfield(pde,'f') && ~isempty(pde.f)
                fq = coefficient(pde.f,xq);
                if size(fq,2) ~= 3
                    error('pde.f must return an N-by-3 body-force array.');
                end
            else
                rho = coefficient(getfielddefault(pde,'rho',1),xq);
                gravity = getfielddefault(pde,'gravity',[0,0,-1]);
                gravity = gravity(:)';
                fq = [rho(:)*gravity(1),rho(:)*gravity(2),rho(:)*gravity(3)];
            end
            for a = 1:10
                F = F+accumarray(double(elem2dof(:,a)),...
                    w(q)*volume.*fq(:,1)*phi(a),[3*Nu,1]);
                F = F+accumarray(Nu+double(elem2dof(:,a)),...
                    w(q)*volume.*fq(:,2)*phi(a),[3*Nu,1]);
                F = F+accumarray(2*Nu+double(elem2dof(:,a)),...
                    w(q)*volume.*fq(:,3)*phi(a),[3*Nu,1]);
            end
        end
    end

    function F = addtraction(F)
        if isempty(topFace) || isempty(pde.g_N), return; end
        [lambdaAll,wAll] = quadpts3face(option.facequadorder);
        faceLocal = ceil((1:size(lambdaAll,1))'/(size(lambdaAll,1)/4));
        for id = 1:numel(topFace)
            [t,f] = ind2sub(size(bdFlag),topFace(id));
            area = localfacearea(t,f);
            localDof = double(elem2dof(t,:));
            for q = find(faceLocal==f)'
                phi = p2value3(lambdaAll(q,:));
                xq = lambdaAll(q,1)*node(elem(t,1),:) + ...
                     lambdaAll(q,2)*node(elem(t,2),:) + ...
                     lambdaAll(q,3)*node(elem(t,3),:) + ...
                     lambdaAll(q,4)*node(elem(t,4),:);
                tq = coefficient(pde.g_N,xq);
                for a = 1:10
                    F(localDof(a)) = F(localDof(a))+wAll(q)*area*tq(1)*phi(a);
                    F(Nu+localDof(a)) = F(Nu+localDof(a))+wAll(q)*area*tq(2)*phi(a);
                    F(2*Nu+localDof(a)) = F(2*Nu+localDof(a))+wAll(q)*area*tq(3)*phi(a);
                end
            end
        end
    end

    function [momentumValue,divergenceValue,constraintValue,...
            totalValue,Kres,Kbres,etaMinRes,etaMaxRes,bedCoefficientRes] = ...
            evaluateresidual(uk,pk,multiplier)
        [Kres,etaMinRes,etaMaxRes] = assembleviscous(uk);
        [Kbres,bedCoefficientRes] = assemblebed(uk,false);
        Fres = assembleforce;
        Fres = addtraction(Fres);
        state = [uk;pk];
        stateMatrix = [Kres+Kbres,B';B,sparse(Np,Np)];
        stateRhs = [Fres;zeros(Np,1)];
        stateResidual = stateMatrix*state-stateRhs+C'*multiplier;
        constraintVector = C*state;
        stateScale = max(1,norm(stateRhs));
        constraintScale = max(1,norm(state));
        momentumValue = norm(stateResidual(1:3*Nu))/stateScale;
        divergenceValue = norm(stateResidual(3*Nu+(1:Np)))/stateScale;
        constraintValue = norm(constraintVector)/constraintScale;
        totalValue = max(norm(stateResidual)/stateScale,constraintValue);
    end

    function [C,nbase] = buildconstraints
        I = [];
        J = [];
        S = [];
        row = 0;
        if option.periodic
            [I,J,S,row] = appendperiodicconstraints(I,J,S,row,udofNode,0);
            [I,J,S,row] = appendperiodicconstraints(I,J,S,row,udofNode,Nu);
            [I,J,S,row] = appendperiodicconstraints(I,J,S,row,udofNode,2*Nu);
            [I,J,S,row] = appendperiodicconstraints(I,J,S,row,node,3*Nu);
        end
        fixedDof = unique(dirichletDof(:));
        if option.periodic
            fixedDof = fixedDof(periodicrepresentativemask(udofNode,fixedDof));
        end
        for c = 0:2
            nr = numel(fixedDof);
            rr = row+(1:nr);
            I = [I,rr]; %#ok<AGROW>
            J = [J,(c*Nu+fixedDof)']; %#ok<AGROW>
            S = [S,ones(1,nr)]; %#ok<AGROW>
            row = row+nr;
        end
        bedOnly = setdiff(unique(bedDof(:)),fixedDof);
        if option.periodic
            bedOnly = bedOnly(periodicrepresentativemask(udofNode,bedOnly));
        end
        nr = numel(bedOnly);
        rr = row+(1:nr);
        normal = bedNormalAtDof(ismember(bedDof,bedOnly),:);
        I = [I,rr,rr,rr];
        J = [J,bedOnly',(Nu+bedOnly)',(2*Nu+bedOnly)'];
        S = [S,normal(:,1)',normal(:,2)',normal(:,3)'];
        row = row+nr;
        nbase = normal;
        if addPressureMeanConstraint
            pressureMean = accumarray(double(elem(:)),...
                repmat(volume/4,4,1),[Np,1]);
            row = row+1;
            I = [I,repmat(row,1,Np)];
            J = [J,3*Nu+(1:Np)];
            S = [S,pressureMean'];
        end
        C = sparse(double(I),double(J),S,row,3*Nu+Np);
    end

    function dofs = boundarydofs(mask)
        dofs = [];
        for id = find(mask(:))'
            [t,f] = ind2sub(size(bdFlag),id);
            localDof = double(elem2dof(t,:));
            activeDof = p2faceactive(f);
            dofs = [dofs;localDof(activeDof)']; %#ok<AGROW>
        end
        dofs = unique(dofs);
    end

    function normal = beddofnormal(dofs)
        normalAccumulator = zeros(numel(dofs),3);
        for id = bedFace'
            [t,f] = ind2sub(size(bdFlag),id);
            faceDof = double(elem2dof(t,p2faceactive(f)));
            faceNormal = localfacenormal(t,f);
            [tf,loc] = ismember(faceDof,dofs);
            normalAccumulator(loc(tf),:) = normalAccumulator(loc(tf),:) + ...
                repmat(faceNormal(:)',sum(tf),1);
        end
        normalNorm = sqrt(sum(normalAccumulator.^2,2));
        normalNorm(normalNorm == 0) = 1;
        normal = normalAccumulator./normalNorm;
    end

    function [I,J,S,row] = appendperiodicconstraints(I,J,S,row,points,offset)
        tol = 1000*eps(max(1,max(abs(points(:)))));
        canonical = periodiccanonical(points,tol);
        key = round(canonical./tol);
        [~,~,group] = unique(key,'rows');
        for ig = 1:max(group)
            members = find(group == ig);
            if numel(members) <= 1
                continue;
            end
            [~,order] = sortrows(points(members,:),[1 2 3]);
            members = members(order);
            representative = members(1);
            slaves = members(2:end);
            for sidx = slaves(:)'
                row = row+1;
                I = [I,row,row]; %#ok<AGROW>
                J = [J,offset+sidx,offset+representative]; %#ok<AGROW>
                S = [S,1,-1]; %#ok<AGROW>
            end
        end
    end

    function mask = periodicrepresentativemask(points,dofs)
        tol = 1000*eps(max(1,max(abs(points(:)))));
        canonical = periodiccanonical(points,tol);
        key = round(canonical./tol);
        [~,~,group] = unique(key,'rows');
        isRepresentative = false(size(points,1),1);
        for ig = 1:max(group)
            members = find(group == ig);
            if isscalar(members)
                isRepresentative(members) = true;
            else
                [~,order] = sortrows(points(members,:),[1 2 3]);
                isRepresentative(members(order(1))) = true;
            end
        end
        mask = isRepresentative(dofs);
    end

    function canonical = periodiccanonical(points,tol)
        canonical = points;
        xl = option.periodic_x(1);
        xr = option.periodic_x(2);
        yl = option.periodic_y(1);
        yr = option.periodic_y(2);
        right = abs(canonical(:,1)-xr) < tol;
        canonical(right,1) = xl;
        topY = abs(canonical(:,2)-yr) < tol;
        canonical(topY,2) = yl;
        canonical(:,3) = canonical(:,3) + ...
            option.periodic_slope_x*points(:,1) + ...
            option.periodic_slope_y*points(:,2);
    end

    function varargout = strainat(ux,uy,uz,Dphi)
        duxdx = sum(Dphi(:,1,:).*reshape(ux(elem2dof),NT,1,10),3);
        duxdy = sum(Dphi(:,2,:).*reshape(ux(elem2dof),NT,1,10),3);
        duxdz = sum(Dphi(:,3,:).*reshape(ux(elem2dof),NT,1,10),3);
        duydx = sum(Dphi(:,1,:).*reshape(uy(elem2dof),NT,1,10),3);
        duydy = sum(Dphi(:,2,:).*reshape(uy(elem2dof),NT,1,10),3);
        duydz = sum(Dphi(:,3,:).*reshape(uy(elem2dof),NT,1,10),3);
        duzdx = sum(Dphi(:,1,:).*reshape(uz(elem2dof),NT,1,10),3);
        duzdy = sum(Dphi(:,2,:).*reshape(uz(elem2dof),NT,1,10),3);
        duzdz = sum(Dphi(:,3,:).*reshape(uz(elem2dof),NT,1,10),3);
        s = {duxdx,duydy,duzdz,...
             0.5*(duxdy+duydx),0.5*(duxdz+duzdx),...
             0.5*(duydz+duzdy)};
        if nargout == 1
            varargout{1} = s;
        else
            varargout = s;
        end
    end

    function area = localfacearea(t,f)
        faceNode = localface(f);
        pts = node(elem(t,faceNode),:);
        area = 0.5*norm(cross(pts(2,:)-pts(1,:),pts(3,:)-pts(1,:)));
    end

    function normal = localfacenormal(t,f)
        faceNode = localface(f);
        pts = node(elem(t,faceNode),:);
        normal = cross(pts(2,:)-pts(1,:),pts(3,:)-pts(1,:));
        if normal(3) < 0
            normal = -normal;
        end
        normal = normal/norm(normal);
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

function phi = p2value3(l)
phi = [l(1)*(2*l(1)-1),l(2)*(2*l(2)-1),...
       l(3)*(2*l(3)-1),l(4)*(2*l(4)-1),...
       4*l(1)*l(2),4*l(1)*l(3),4*l(1)*l(4),...
       4*l(2)*l(3),4*l(2)*l(4),4*l(3)*l(4)];
end

function xq = barycenter(l,node,elem)
xq = l(1)*node(elem(:,1),:) + l(2)*node(elem(:,2),:) + ...
     l(3)*node(elem(:,3),:) + l(4)*node(elem(:,4),:);
end

function [component,basis] = splitlocal(a)
component = ceil(a/10);
basis = a-10*(component-1);
end

function s = basisstrain3(component,Dphi,basis)
z = zeros(size(Dphi,1),1);
s = {z,z,z,z,z,z};
if component == 1
    s{1} = Dphi(:,1,basis);
    s{4} = 0.5*Dphi(:,2,basis);
    s{5} = 0.5*Dphi(:,3,basis);
elseif component == 2
    s{2} = Dphi(:,2,basis);
    s{4} = 0.5*Dphi(:,1,basis);
    s{6} = 0.5*Dphi(:,3,basis);
else
    s{3} = Dphi(:,3,basis);
    s{5} = 0.5*Dphi(:,1,basis);
    s{6} = 0.5*Dphi(:,2,basis);
end
end

function value = straininner(a,b)
value = a{1}.*b{1}+a{2}.*b{2}+a{3}.*b{3}+...
    2*(a{4}.*b{4}+a{5}.*b{5}+a{6}.*b{6});
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

function value = getfielddefault(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end

function face = localface(f)
faces = [2 4 3;1 3 4;1 4 2;1 2 3];
face = faces(f,:);
end

function active = p2faceactive(zeroLambda)
active = true(1,10);
if zeroLambda == 1
    active([1 5 6 7]) = false;
elseif zeroLambda == 2
    active([2 5 8 9]) = false;
elseif zeroLambda == 3
    active([3 6 8 10]) = false;
else
    active([4 7 9 10]) = false;
end
end
