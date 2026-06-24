function [soln,eqn,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option)
%% NONLINEARSTOKESP2P1_PERIODIC Nonlinear full-Stokes ice-flow model.
%
% Solve, in a two-dimensional x-z cross-section,
%
%   div(u) = 0,
%   -div(eta(u)*(grad(u)+grad(u)')) + grad(p) = f,
%
% with P2 velocity and P1 pressure.  Boundary flags are
%
%   2: prescribed traction sigma*n = g_N (zero when pde.g_N is absent),
%   3: impermeable Weertman bed
%          u*n = 0,
%          T*sigma*n + beta*|T*u|^(m-1)*T*u = 0.
%
% If option.periodic is true, x = option.periodic_x(1) and
% x = option.periodic_x(2) are identified.  The two boundaries may be
% vertically shifted, as in a parallelogram ice slab.
%
% Glen's law is
%
%   eta = 0.5*A^(-1/n)*(epsII + eps_reg^2)^((1-n)/(2*n)),
%   epsII = 0.5*trace(eps(u)^2),
%   eps(u) = 0.5*(grad(u)+grad(u)').
%
% Required data:
%   pde.A, pde.n, pde.beta, pde.m
%
% Optional data:
%   pde.f, pde.rho, pde.gravity, pde.g_N
%
% Important options:
%   option.periodic, option.periodic_x, option.tol,
%   option.maxIt, option.damping, option.eps_reg, option.quadorder,
%   option.pressure_constraint ('auto', 'mean-zero', or 'none'),
%   option.residual_tol, option.residual_check_threshold.
%
% The nonlinear iteration is a damped Picard iteration.  The small
% eps_reg regularizes Glen's singular viscosity at zero strain rate.

if nargin < 5, option = struct; end
option = setoption(option,'periodic',true);
option = setoption(option,'periodic_x',[min(node(:,1)),max(node(:,1))]);
option = setoption(option,'tol',1e-8);
option = setoption(option,'maxIt',50);
option = setoption(option,'damping',0.7);
option = setoption(option,'eps_reg',1e-8);
option = setoption(option,'quadorder',4);
option = setoption(option,'printlevel',1);
option = setoption(option,'assemble_tangent',false);
option = setoption(option,'pressure_constraint','auto');
option = setoption(option,'residual_tol',option.tol);
option = setoption(option,'residual_check_threshold',...
    max(1e-2,sqrt(option.residual_tol)));

if ~isfield(pde,'A'), pde.A = 1; end
if ~isfield(pde,'n'), pde.n = 3; end
if ~isfield(pde,'beta'), pde.beta = 1; end
if ~isfield(pde,'m'), pde.m = 1; end
if ~isfield(pde,'g_N'), pde.g_N = []; end

[elem2dof,edge] = dofP2(elem);
[Dlambda,area] = gradbasis(node,elem);
N = size(node,1);
NT = size(elem,1);
NE = size(edge,1);
Nu = N + NE;
Np = N;
udofNode = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];

elem2edge = elem2dof(:,4:6)-N;
isBedEdge = false(NE,1);
isTopEdge = false(NE,1);
if ~isempty(bdFlag)
    isBedEdge(elem2edge(bdFlag(:)==3)) = true;
    isTopEdge(elem2edge(bdFlag(:)==2)) = true;
end
bedEdgeIdx = find(isBedEdge);
topEdgeIdx = find(isTopEdge);

hasTractionBoundary = any(bdFlag(:)==2);
if ~(ischar(option.pressure_constraint) || ...
        (isstring(option.pressure_constraint) && isscalar(option.pressure_constraint)))
    error('iFEM:NonlinearStokesPressureConstraint',...
        'pressure_constraint must be a character vector or scalar string.');
end
pressureConstraint = lower(strtrim(char(option.pressure_constraint)));
switch pressureConstraint
    case 'auto'
        addPressureMeanConstraint = ~hasTractionBoundary;
    case 'mean-zero'
        addPressureMeanConstraint = true;
    case 'none'
        addPressureMeanConstraint = false;
    otherwise
        error('iFEM:NonlinearStokesPressureConstraint',...
            'Unknown pressure_constraint value: %s.',pressureConstraint);
end
if ~(isscalar(option.residual_tol) && isfinite(option.residual_tol) && ...
        option.residual_tol > 0)
    error('iFEM:NonlinearStokesResidualTolerance',...
        'residual_tol must be a positive finite scalar.');
end
if ~(isscalar(option.residual_check_threshold) && ...
        isfinite(option.residual_check_threshold) && ...
        option.residual_check_threshold > 0)
    error('iFEM:NonlinearStokesResidualThreshold',...
        'residual_check_threshold must be a positive finite scalar.');
end

[C,bedDof,bedNormal] = buildconstraints;
nConstraint = size(C,1);
B = Bmatrix;

u = zeros(2*Nu,1);
p = zeros(Np,1);
constraintMultiplier = zeros(nConstraint,1);
if isfield(option,'u0') && numel(option.u0) == 2*Nu
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
    [Kb,bedCoefficient] = assemblebed(u);
    F = assembleforce;
    F = addtraction(F);

    M = [K+Kb, B'; B, sparse(Np,Np)];
    rhs = [F; zeros(Np,1)];
    saddle = [M, C'; C, sparse(nConstraint,nConstraint)];
    fullsol = saddle\[rhs; zeros(nConstraint,1)];

    unew = fullsol(1:2*Nu);
    pnew = fullsol(2*Nu+(1:Np));
    multiplierNew = fullsol(2*Nu+Np+(1:nConstraint));
    alpha = option.damping;
    updatedU = (1-alpha)*u + alpha*unew;
    updatedP = (1-alpha)*p + alpha*pnew;
    updatedMultiplier = (1-alpha)*constraintMultiplier + ...
        alpha*multiplierNew;

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
            fprintf(['nonlinear Stokes %2d: relchange %.3e, ',...
                'residual %.3e, eta [%.3e, %.3e]\n'],...
                k,residual(k),nonlinearResidual(k),etaMin,etaMax);
        else
            fprintf(['nonlinear Stokes %2d: relchange %.3e, ',...
                'eta [%.3e, %.3e]\n'],...
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
info.hasTractionBoundary = hasTractionBoundary;
info.pressureMeanConstrained = addPressureMeanConstraint;

soln.u = u;
soln.p = p;
soln.ux = u(1:Nu);
soln.uz = u(Nu+1:end);

eqn.K = K+Kb;
eqn.B = B;
eqn.C = C;
eqn.edge = edge;
eqn.elem2dof = elem2dof;
eqn.bedDof = bedDof;
eqn.bedNormal = bedNormal;
eqn.constraintMultiplier = constraintMultiplier;

if option.assemble_tangent
    % Consistent Jacobian of the converged nonlinear residual.  It is
    % required by incremental-state and adjoint equations; the Picard
    % matrix K+Kb omits derivatives of viscosity and sliding coefficient.
    Kt = assembleviscoustangent(u);
    Kbt = assemblebedtangent(u);
    tangentM = [Kt+Kbt, B'; B, sparse(Np,Np)];
    eqn.tangent = [tangentM, C'; C, sparse(nConstraint,nConstraint)];
    eqn.applyBetaDerivative = @assemblebetadirection;
end

    function [K,etaMin,etaMax] = assembleviscous(uk)
        [lambda,w] = quadpts(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*144*NT,1);
        cols = zeros(nq*144*NT,1);
        vals = zeros(nq*144*NT,1);
        cursor = 0;
        etaMin = inf;
        etaMax = 0;

        ux = uk(1:Nu);
        uz = uk(Nu+1:end);
        for q = 1:nq
            Dphi = p2gradient(lambda(q,:),Dlambda);
            duxdx = sum(Dphi(:,1,:).*reshape(ux(elem2dof),NT,1,6),3);
            duxdz = sum(Dphi(:,2,:).*reshape(ux(elem2dof),NT,1,6),3);
            duzdx = sum(Dphi(:,1,:).*reshape(uz(elem2dof),NT,1,6),3);
            duzdz = sum(Dphi(:,2,:).*reshape(uz(elem2dof),NT,1,6),3);
            exx = duxdx;
            ezz = duzdz;
            exz = 0.5*(duxdz+duzdx);
            epsII = 0.5*(exx.^2+ezz.^2+2*exz.^2);
            xq = lambda(q,1)*node(elem(:,1),:) + ...
                 lambda(q,2)*node(elem(:,2),:) + ...
                 lambda(q,3)*node(elem(:,3),:);
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            eta = 0.5.*Aq.^(-1./nqfield).*...
                (epsII+option.eps_reg^2).^((1-nqfield)./(2*nqfield));
            etaMin = min(etaMin,min(eta));
            etaMax = max(etaMax,max(eta));

            localDof = [elem2dof, Nu+elem2dof];
            for a = 1:12
                [aComp,aBasis] = splitlocal(a);
                [aExx,aEzz,aExz] = basisstrain(aComp,Dphi,aBasis);
                for b = 1:12
                    [bComp,bBasis] = splitlocal(b);
                    [bExx,bEzz,bExz] = basisstrain(bComp,Dphi,bBasis);
                    kab = 2*w(q)*area.*eta.*...
                        (aExx.*bExx+aEzz.*bEzz+2*aExz.*bExz);
                    idx = cursor+(1:NT);
                    rows(idx) = localDof(:,a);
                    cols(idx) = localDof(:,b);
                    vals(idx) = kab;
                    cursor = cursor+NT;
                end
            end
        end
        K = sparse(rows,cols,vals,2*Nu,2*Nu);
    end

    function B = Bmatrix
        [lambda,w] = quadpts(2);
        Bx = sparse(Np,Nu);
        Bz = sparse(Np,Nu);
        for q = 1:size(lambda,1)
            Dphi = p2gradient(lambda(q,:),Dlambda);
            for a = 1:6
                for j = 1:3
                    v = w(q)*area.*lambda(q,j);
                    Bx = Bx+sparse(double(elem(:,j)),double(elem2dof(:,a)),...
                        -v.*Dphi(:,1,a),Np,Nu);
                    Bz = Bz+sparse(double(elem(:,j)),double(elem2dof(:,a)),...
                        -v.*Dphi(:,2,a),Np,Nu);
                end
            end
        end
        B = [Bx,Bz];
    end

    function [Kb,coefAtMid] = assemblebed(uk)
        Kb = sparse(2*Nu,2*Nu);
        coefAtMid = [];
        if isempty(bedEdgeIdx), return; end
        [lbd,wbd] = quadpts1(6);
        bed = edge(bedEdgeIdx,:);
        bedLocalDof = [bed, N+bedEdgeIdx];
        tangent = node(bed(:,2),:)-node(bed(:,1),:);
        edgeLength = sqrt(sum(tangent.^2,2));
        tangent = tangent./edgeLength;
        ubx = uk(bedLocalDof);
        ubz = uk(Nu+bedLocalDof);
        ii = [];
        jj = [];
        ss = [];
        coefAtMid = zeros(length(bedEdgeIdx),1);
        for q = 1:size(lbd,1)
            phi = [(2*lbd(q,1)-1)*lbd(q,1),...
                   (2*lbd(q,2)-1)*lbd(q,2),...
                   4*lbd(q,1)*lbd(q,2)];
            xq = lbd(q,1)*node(bed(:,1),:)+lbd(q,2)*node(bed(:,2),:);
            beta = coefficient(pde.beta,xq);
            ut = sum((ubx*phi').*tangent(:,1)+(ubz*phi').*tangent(:,2),2);
            gamma = beta.*(ut.^2+option.eps_reg^2).^((pde.m-1)/2);
            if q == ceil(size(lbd,1)/2), coefAtMid = gamma; end
            for a = 1:3
                ia = bedLocalDof(:,a);
                for b = 1:3
                    ib = bedLocalDof(:,b);
                    s = wbd(q)*edgeLength.*gamma*phi(a)*phi(b);
                    ii = [ii;ia;ia;Nu+ia;Nu+ia]; %#ok<AGROW>
                    jj = [jj;ib;Nu+ib;ib;Nu+ib]; %#ok<AGROW>
                    ss = [ss;s.*tangent(:,1).^2;...
                        s.*tangent(:,1).*tangent(:,2);...
                        s.*tangent(:,2).*tangent(:,1);...
                        s.*tangent(:,2).^2]; %#ok<AGROW>
                end
            end
        end
        Kb = sparse(double(ii),double(jj),ss,2*Nu,2*Nu);
    end

    function Kt = assembleviscoustangent(uk)
        [lambda,w] = quadpts(option.quadorder);
        nq = size(lambda,1);
        rows = zeros(nq*144*NT,1);
        cols = zeros(nq*144*NT,1);
        vals = zeros(nq*144*NT,1);
        cursor = 0;
        ux = uk(1:Nu);
        uz = uk(Nu+1:end);

        for q = 1:nq
            Dphi = p2gradient(lambda(q,:),Dlambda);
            duxdx = sum(Dphi(:,1,:).*reshape(ux(elem2dof),NT,1,6),3);
            duxdz = sum(Dphi(:,2,:).*reshape(ux(elem2dof),NT,1,6),3);
            duzdx = sum(Dphi(:,1,:).*reshape(uz(elem2dof),NT,1,6),3);
            duzdz = sum(Dphi(:,2,:).*reshape(uz(elem2dof),NT,1,6),3);
            exx = duxdx;
            ezz = duzdz;
            exz = 0.5*(duxdz+duzdx);
            epsII = 0.5*(exx.^2+ezz.^2+2*exz.^2);
            xq = lambda(q,1)*node(elem(:,1),:) + ...
                 lambda(q,2)*node(elem(:,2),:) + ...
                 lambda(q,3)*node(elem(:,3),:);
            Aq = coefficient(pde.A,xq);
            nqfield = coefficient(pde.n,xq);
            strainRegularized = epsII+option.eps_reg^2;
            exponent = (1-nqfield)./(2*nqfield);
            eta = 0.5.*Aq.^(-1./nqfield).*...
                strainRegularized.^exponent;

            localDof = [elem2dof,Nu+elem2dof];
            for a = 1:12
                [aComp,aBasis] = splitlocal(a);
                [aExx,aEzz,aExz] = basisstrain(aComp,Dphi,aBasis);
                stateDotA = exx.*aExx+ezz.*aEzz+2*exz.*aExz;
                for b = 1:12
                    [bComp,bBasis] = splitlocal(b);
                    [bExx,bEzz,bExz] = basisstrain(bComp,Dphi,bBasis);
                    stateDotB = exx.*bExx+ezz.*bEzz+2*exz.*bExz;
                    strainDot = aExx.*bExx+aEzz.*bEzz+...
                        2*aExz.*bExz;
                    kab = 2*w(q)*area.*eta.*...
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
        Kt = sparse(rows,cols,vals,2*Nu,2*Nu);
    end

    function Kbt = assemblebedtangent(uk)
        Kbt = sparse(2*Nu,2*Nu);
        if isempty(bedEdgeIdx), return; end
        [lbd,wbd] = quadpts1(6);
        bed = edge(bedEdgeIdx,:);
        bedLocalDof = [bed,N+bedEdgeIdx];
        tangent = node(bed(:,2),:)-node(bed(:,1),:);
        edgeLength = sqrt(sum(tangent.^2,2));
        tangent = tangent./edgeLength;
        ubx = uk(bedLocalDof);
        ubz = uk(Nu+bedLocalDof);
        ii = [];
        jj = [];
        ss = [];
        exponent = (pde.m-1)/2;

        for q = 1:size(lbd,1)
            phi = [(2*lbd(q,1)-1)*lbd(q,1),...
                   (2*lbd(q,2)-1)*lbd(q,2),...
                   4*lbd(q,1)*lbd(q,2)];
            xq = lbd(q,1)*node(bed(:,1),:)+...
                 lbd(q,2)*node(bed(:,2),:);
            beta = coefficient(pde.beta,xq);
            ut = sum((ubx*phi').*tangent(:,1)+...
                     (ubz*phi').*tangent(:,2),2);
            speedRegularized = ut.^2+option.eps_reg^2;
            tangentCoefficient = beta.*...
                (speedRegularized.^exponent+...
                 (pde.m-1)*ut.^2.*speedRegularized.^(exponent-1));
            for a = 1:3
                ia = bedLocalDof(:,a);
                for b = 1:3
                    ib = bedLocalDof(:,b);
                    s = wbd(q)*edgeLength.*tangentCoefficient*...
                        phi(a)*phi(b);
                    ii = [ii;ia;ia;Nu+ia;Nu+ia]; %#ok<AGROW>
                    jj = [jj;ib;Nu+ib;ib;Nu+ib]; %#ok<AGROW>
                    ss = [ss;s.*tangent(:,1).^2;...
                        s.*tangent(:,1).*tangent(:,2);...
                        s.*tangent(:,2).*tangent(:,1);...
                        s.*tangent(:,2).^2]; %#ok<AGROW>
                end
            end
        end
        Kbt = sparse(double(ii),double(jj),ss,2*Nu,2*Nu);
    end

    function load = assemblebetadirection(betaDirection)
        % Derivative of the nonlinear residual for a supplied delta-beta.
        load = zeros(2*Nu+Np+nConstraint,1);
        if isempty(bedEdgeIdx), return; end
        [lbd,wbd] = quadpts1(6);
        bed = edge(bedEdgeIdx,:);
        bedLocalDof = [bed,N+bedEdgeIdx];
        tangent = node(bed(:,2),:)-node(bed(:,1),:);
        edgeLength = sqrt(sum(tangent.^2,2));
        tangent = tangent./edgeLength;
        ubx = u(bedLocalDof);
        ubz = u(Nu+bedLocalDof);
        exponent = (pde.m-1)/2;

        for q = 1:size(lbd,1)
            phi = [(2*lbd(q,1)-1)*lbd(q,1),...
                   (2*lbd(q,2)-1)*lbd(q,2),...
                   4*lbd(q,1)*lbd(q,2)];
            xq = lbd(q,1)*node(bed(:,1),:)+...
                 lbd(q,2)*node(bed(:,2),:);
            deltaBeta = coefficient(betaDirection,xq);
            ut = sum((ubx*phi').*tangent(:,1)+...
                     (ubz*phi').*tangent(:,2),2);
            tractionDirection = deltaBeta.*...
                (ut.^2+option.eps_reg^2).^exponent.*ut;
            for a = 1:3
                contribution = wbd(q)*edgeLength.*...
                    tractionDirection*phi(a);
                ia = bedLocalDof(:,a);
                load(1:2*Nu) = load(1:2*Nu)+...
                    accumarray(ia,contribution.*tangent(:,1),...
                               [2*Nu,1]);
                load(1:2*Nu) = load(1:2*Nu)+...
                    accumarray(Nu+ia,contribution.*tangent(:,2),...
                               [2*Nu,1]);
            end
        end
    end

    function F = assembleforce
        F = zeros(2*Nu,1);
        [lambda,w] = quadpts(option.quadorder);
        for q = 1:size(lambda,1)
            phi = [lambda(q,1)*(2*lambda(q,1)-1),...
                   lambda(q,2)*(2*lambda(q,2)-1),...
                   lambda(q,3)*(2*lambda(q,3)-1),...
                   4*lambda(q,2)*lambda(q,3),...
                   4*lambda(q,3)*lambda(q,1),...
                   4*lambda(q,1)*lambda(q,2)];
            xq = lambda(q,1)*node(elem(:,1),:) + ...
                 lambda(q,2)*node(elem(:,2),:) + ...
                 lambda(q,3)*node(elem(:,3),:);
            if isfield(pde,'f') && ~isempty(pde.f)
                fq = coefficient(pde.f,xq);
                if size(fq,2) ~= 2
                    error('pde.f must return an N-by-2 body-force array.');
                end
            else
                rho = coefficient(getfielddefault(pde,'rho',1),xq);
                gravity = getfielddefault(pde,'gravity',[0,-1]);
                if numel(gravity) ~= 2
                    error('pde.gravity must contain the x and z components.');
                end
                gravity = gravity(:)';
                fq = [rho(:)*gravity(1),rho(:)*gravity(2)];
            end
            for a = 1:6
                F = F+accumarray(elem2dof(:,a),...
                    w(q)*area.*fq(:,1)*phi(a),[2*Nu,1]);
                F = F+accumarray(Nu+elem2dof(:,a),...
                    w(q)*area.*fq(:,2)*phi(a),[2*Nu,1]);
            end
        end
    end

    function [momentumValue,divergenceValue,constraintValue,...
            totalValue,Kres,Kbres,etaMinRes,etaMaxRes,bedCoefficientRes] = ...
            evaluateresidual(uk,pk,multiplier)
        [Kres,etaMinRes,etaMaxRes] = assembleviscous(uk);
        [Kbres,bedCoefficientRes] = assemblebed(uk);
        Fres = assembleforce;
        Fres = addtraction(Fres);
        state = [uk;pk];
        stateMatrix = [Kres+Kbres,B';B,sparse(Np,Np)];
        stateRhs = [Fres;zeros(Np,1)];
        stateResidual = stateMatrix*state-stateRhs+C'*multiplier;
        constraintVector = C*state;
        stateScale = max(1,norm(stateRhs));
        constraintScale = max(1,norm(state));
        momentumValue = norm(stateResidual(1:2*Nu))/stateScale;
        divergenceValue = norm(stateResidual(2*Nu+(1:Np)))/stateScale;
        constraintValue = norm(constraintVector)/constraintScale;
        totalValue = max(norm(stateResidual)/stateScale,constraintValue);
    end

    function F = addtraction(F)
        if isempty(topEdgeIdx) || isempty(pde.g_N), return; end
        [lbd,wbd] = quadpts1(6);
        top = edge(topEdgeIdx,:);
        topDof = [top,N+topEdgeIdx];
        edgeLength = sqrt(sum((node(top(:,2),:)-node(top(:,1),:)).^2,2));
        for q = 1:size(lbd,1)
            phi = [(2*lbd(q,1)-1)*lbd(q,1),...
                   (2*lbd(q,2)-1)*lbd(q,2),...
                   4*lbd(q,1)*lbd(q,2)];
            xq = lbd(q,1)*node(top(:,1),:)+lbd(q,2)*node(top(:,2),:);
            tq = coefficient(pde.g_N,xq);
            for a = 1:3
                F = F+accumarray(topDof(:,a),...
                    wbd(q)*edgeLength.*tq(:,1)*phi(a),[2*Nu,1]);
                F = F+accumarray(Nu+topDof(:,a),...
                    wbd(q)*edgeLength.*tq(:,2)*phi(a),[2*Nu,1]);
            end
        end
    end

    function [C,baseDof,nbase] = buildconstraints
        I = [];
        J = [];
        S = [];
        row = 0;
        tolx = 100*eps(max(1,max(abs(node(:)))));
        if option.periodic
            xl = option.periodic_x(1);
            xr = option.periodic_x(2);
            leftU = find(abs(udofNode(:,1)-xl)<tolx);
            rightU = find(abs(udofNode(:,1)-xr)<tolx);
            leftU = sortperiodic(leftU,udofNode);
            rightU = sortperiodic(rightU,udofNode);
            assert(length(leftU)==length(rightU),'Periodic velocity boundaries do not match.');
            for c = 0:1
                nr = length(leftU);
                rr = row+(1:nr);
                I = [I,rr,rr]; %#ok<AGROW>
                J = [J,(leftU+c*Nu)',(rightU+c*Nu)']; %#ok<AGROW>
                S = [S,ones(1,nr),-ones(1,nr)]; %#ok<AGROW>
                row = row+nr;
            end
            leftP = find(abs(node(:,1)-xl)<tolx);
            rightP = find(abs(node(:,1)-xr)<tolx);
            leftP = sortperiodic(leftP,node);
            rightP = sortperiodic(rightP,node);
            assert(length(leftP)==length(rightP),'Periodic pressure boundaries do not match.');
            nr = length(leftP);
            rr = row+(1:nr);
            I = [I,rr,rr];
            J = [J,(2*Nu+leftP)',(2*Nu+rightP)'];
            S = [S,ones(1,nr),-ones(1,nr)];
            row = row+nr;
        end

        baseDof = unique([edge(bedEdgeIdx,:)';(N+bedEdgeIdx)']);
        xb = udofNode(baseDof,:);
        [~,order] = sort(xb(:,1));
        baseDof = baseDof(order);
        if option.periodic
            baseDof(abs(udofNode(baseDof,1)-option.periodic_x(2))<tolx) = [];
        end
        nbase = zeros(length(baseDof),2);
        for ib = 1:length(baseDof)
            attached = bedEdgeIdx(any(edge(bedEdgeIdx,:)==baseDof(ib),2));
            if baseDof(ib)>N
                attached = baseDof(ib)-N;
            end
            tangent = node(edge(attached,2),:)-node(edge(attached,1),:);
            tangent = sum(tangent./sqrt(sum(tangent.^2,2)),1);
            tangent = tangent/norm(tangent);
            nbase(ib,:) = [tangent(2),-tangent(1)];
        end
        nr = length(baseDof);
        rr = row+(1:nr);
        I = [I,rr,rr];
        J = [J,baseDof',(Nu+baseDof)'];
        S = [S,nbase(:,1)',nbase(:,2)'];
        row = row+nr;

        if addPressureMeanConstraint
            pressureMean = accumarray(double(elem(:)),...
                repmat(area/3,3,1),[Np,1]);
            row = row+1;
            I = [I,repmat(row,1,Np)];
            J = [J,2*Nu+(1:Np)];
            S = [S,pressureMean'];
        end
        C = sparse(double(I),double(J),S,row,2*Nu+Np);
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

function [component,basis] = splitlocal(a)
component = 1+(a>6);
basis = a-6*(component-1);
end

function [exx,ezz,exz] = basisstrain(component,Dphi,basis)
z = zeros(size(Dphi,1),1);
if component == 1
    exx = Dphi(:,1,basis);
    ezz = z;
    exz = 0.5*Dphi(:,2,basis);
else
    exx = z;
    ezz = Dphi(:,2,basis);
    exz = 0.5*Dphi(:,1,basis);
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

function value = getfielddefault(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end

function idx = sortperiodic(idx,points)
[~,order] = sort(points(idx,2));
idx = idx(order);
end
