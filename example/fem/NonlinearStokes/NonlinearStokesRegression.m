%% NONLINEARSTOKESREGRESSION Full nonlinear Stokes regression suite.
%
% Run from the iFEM root after setpath:
%
%   NonlinearStokesRegression
%
% The suite checks the forward solver, pressure constraints, manufactured
% convergence rates, continuation, adjoint derivatives, and inversion
% stability.  It is intentionally a multi-minute test.

close all;
clear variables;
set(groot,'defaultFigureVisible','off');

fprintf('=== nonlinear Stokes regression ===\n');

%% Ice slab: grid convergence, mass conservation, and full residual
hlist = [1/8;1/16;1/32];
divergenceL2 = zeros(size(hlist));
ice = cell(size(hlist));
for level = 1:numel(hlist)
    [ice{level},node,elem] = solveiceslab(...
        hlist(level),'auto',1e-4,[],1e-8);
    divergenceL2(level) = getdivergenceL2(...
        node,elem,ice{level}.soln);
    assert(ice{level}.info.converged,...
        'Ice-slab solve failed on h=%g.',hlist(level));
    assert(~ice{level}.info.pressureMeanConstrained,...
        'Auto mode incorrectly constrained pressure with traction data.');
    assert(ice{level}.info.nonlinearResidual(end) < 1e-8,...
        'Ice-slab nonlinear residual is above tolerance.');
    assert(abs(getdivergenceintegral(node,elem,ice{level}.soln)) < 1e-12,...
        'Ice-slab global mass balance failed.');
    fprintf('ice h=%7.5f: divL2 %.3e, residual %.3e, it %d\n',...
        hlist(level),divergenceL2(level),...
        ice{level}.info.nonlinearResidual(end),ice{level}.info.itStep);
end
assert(all(diff(divergenceL2)<0),...
    'Ice-slab divergence did not decrease under mesh refinement.');

%% Pressure-constraint modes
[iceNone,~,~] = solveiceslab(1/8,'none',1e-4,[],1e-8);
[iceMean,nodeMean,elemMean] = solveiceslab(...
    1/8,'mean-zero',1e-4,[],1e-8);
assert(norm(ice{1}.soln.u-iceNone.soln.u)/norm(ice{1}.soln.u) < 1e-10,...
    'Auto and none velocity solutions differ with traction data.');
assert(norm(ice{1}.soln.p-iceNone.soln.p)/norm(ice{1}.soln.p) < 1e-10,...
    'Auto and none pressure solutions differ with traction data.');
assert(iceMean.info.pressureMeanConstrained,...
    'Explicit mean-zero mode did not add the pressure constraint.');
assert(getdivergenceL2(nodeMean,elemMean,iceMean.soln) > 1e-3,...
    'Mean-zero ice-slab diagnostic no longer exposes incompatibility.');

invalidModeRejected = false;
try
    solveiceslab(1/4,'invalid-mode',1e-4,[],1e-8);
catch exception
    invalidModeRejected = strcmp(exception.identifier,...
        'iFEM:NonlinearStokesPressureConstraint');
end
assert(invalidModeRejected,'Unknown pressure-constraint mode was accepted.');

%% Manufactured-solution convergence
L = 1;
H = 1;
slope = 0.1;
epsReg = 1e-2;
hMms = [1/4;1/8;1/16];
errU = zeros(size(hMms));
errP = zeros(size(hMms));
divMms = zeros(size(hMms));
for level = 1:numel(hMms)
    pde = NonlinearStokesMMSData(L,H,slope,epsReg);
    [node,elem] = squaremesh([0,L,0,H],hMms(level));
    bdFlag = setboundary(node,elem,...
        'Neumann','y==1','Robin','y==0');
    node(:,2) = node(:,2)-slope*node(:,1);
    option = forwardoption(L,epsReg,1e-11);
    option.maxIt = 180;
    option.quadorder = 9;
    [soln,~,info] = NonlinearStokesP2P1(...
        node,elem,bdFlag,pde,option);
    assert(info.converged,'MMS solve failed on h=%g.',hMms(level));
    errUx = getL2error(node,elem,pde.exactux,soln.ux,6);
    errUz = getL2error(node,elem,pde.exactuz,soln.uz,6);
    errU(level) = hypot(errUx,errUz);
    errP(level) = getL2error(node,elem,pde.exactp,soln.p,6);
    divMms(level) = getdivergenceL2(node,elem,soln);
end
rateU = log(errU(1:end-1)./errU(2:end))/log(2);
rateP = log(errP(1:end-1)./errP(2:end))/log(2);
assert(rateU(end)>2.7,'MMS velocity convergence rate is too low.');
assert(rateP(end)>1.7,'MMS pressure convergence rate is too low.');
assert(all(diff(divMms)<0),'MMS divergence did not decrease.');
fprintf('MMS finest rates: velocity %.3f, pressure %.3f\n',...
    rateU(end),rateP(end));

% The manufactured pressure already has zero mean, so both modes agree.
pde = NonlinearStokesMMSData(L,H,slope,epsReg);
[node,elem] = squaremesh([0,L,0,H],1/8);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);
option = forwardoption(L,epsReg,1e-11);
option.maxIt = 180;
option.quadorder = 9;
[mmsAuto,~,autoInfo] = NonlinearStokesP2P1(...
    node,elem,bdFlag,pde,option);
option.pressure_constraint = 'mean-zero';
[mmsMean,~,meanInfo] = NonlinearStokesP2P1(...
    node,elem,bdFlag,pde,option);
assert(~autoInfo.pressureMeanConstrained && ...
       meanInfo.pressureMeanConstrained,...
    'MMS pressure modes were not selected as requested.');
assert(norm(mmsAuto.u-mmsMean.u)/norm(mmsAuto.u)<1e-10,...
    'MMS velocity changed under compatible pressure normalization.');
assert(norm(mmsAuto.p-mmsMean.p)/norm(mmsAuto.p)<1e-10,...
    'MMS pressure changed under compatible pressure normalization.');

%% No-traction case: auto mode must remove the pressure nullspace
[node,elem] = squaremesh([0,1,0,1],1/4);
bdFlag = setboundary(node,elem,'Robin','all');
pde = struct('A',1,'n',1,'beta',1,'m',1,...
    'rho',1,'gravity',[0,0],'g_N',[]);
option = forwardoption(1,1e-4,1e-10);
option.periodic = false;
option.maxIt = 20;
[soln,eqn,info] = NonlinearStokesP2P1(...
    node,elem,bdFlag,pde,option);
assert(info.converged,'No-traction pressure-nullspace test did not converge.');
assert(~info.hasTractionBoundary && info.pressureMeanConstrained,...
    'Auto mode did not constrain pressure in the no-traction case.');
pressureMass = getpressuremass(node,elem);
assert(abs(pressureMass'*soln.p)<1e-12,...
    'No-traction pressure does not satisfy the zero-mean constraint.');
Np = numel(soln.p);
stateMatrix = [eqn.K,eqn.B';eqn.B,sparse(Np,Np)];
saddle = [stateMatrix,eqn.C';...
          eqn.C,sparse(size(eqn.C,1),size(eqn.C,1))];
assert(sprank(saddle)==size(saddle,1),...
    'No-traction augmented matrix is structurally singular.');

%% Direct solve versus regularization continuation
warm = [];
for epsValue = 10.^(-1:-1:-4)
    [stage,~,~] = solveiceslab(1/8,'auto',epsValue,warm,1e-8);
    assert(stage.info.converged,...
        'Continuation failed at eps_reg=%g.',epsValue);
    warm = stage.soln.u;
end
assert(norm(stage.soln.u-ice{1}.soln.u)/norm(ice{1}.soln.u)<2e-6,...
    'Direct and continuation velocities differ.');
assert(norm(stage.soln.p-ice{1}.soln.p)/norm(ice{1}.soln.p)<1e-8,...
    'Direct and continuation pressures differ.');

fprintf('Forward, pressure, MMS, and continuation checks passed.\n');

%% Adjoint inversion and derivative checks
NonlinearStokesAdjointInversion;
assert(all(isfinite(history.objective)),...
    'Adjoint inversion objective contains NaN or Inf.');
assert(all(diff(history.objective)<=1e-12),...
    'Adjoint inversion objective is not monotonically decreasing.');
assert(derivativeCheck.stateError<1e-4,...
    'Incremental-state derivative check failed.');
assert(derivativeCheck.gradientError<1e-4,...
    'Adjoint-gradient derivative check failed.');
assert(derivativeCheck.gaussNewtonError<1e-4,...
    'Gauss-Newton derivative check failed.');

%% Finite-difference inversion stability
NonlinearStokesBetaInversion;
assert(all(isfinite(history.objective)),...
    'Finite-difference inversion objective contains NaN or Inf.');
assert(all(isfinite(history.dataMisfit)),...
    'Finite-difference inversion data misfit contains NaN or Inf.');
assert(all(diff(history.objective)<=1e-12),...
    'Finite-difference inversion objective is not monotonically decreasing.');
assert(all(isfinite(betaRecovered)) && all(betaRecovered>0),...
    'Finite-difference inversion produced an invalid beta.');

fprintf('All nonlinear Stokes regression checks passed.\n');

function option = forwardoption(L,epsReg,tolerance)
    option.periodic = true;
    option.periodic_x = [0,L];
    option.eps_reg = epsReg;
    option.maxIt = 150;
    option.tol = tolerance;
    option.residual_tol = tolerance;
    option.damping = 0.8;
    option.printlevel = 0;
    option.quadorder = 6;
end

function [result,node,elem] = solveiceslab(...
        h,pressureMode,epsReg,u0,tolerance)
    L = 1;
    H = 1;
    slope = 0.1;
    [node,elem] = squaremesh([0,L,0,H],h);
    bdFlag = setboundary(node,elem,...
        'Neumann','y==1','Robin','y==0');
    node(:,2) = node(:,2)-slope*node(:,1);
    pde = struct('A',1,'n',3,'beta',10,'m',1/3,...
        'rho',1,'gravity',[0,-1],'g_N',[]);
    option = forwardoption(L,epsReg,tolerance);
    option.maxIt = 180;
    option.pressure_constraint = pressureMode;
    if ~isempty(u0)
        option.u0 = u0;
    end
    [result.soln,result.eqn,result.info] = ...
        NonlinearStokesP2P1(node,elem,bdFlag,pde,option);
end

function value = getdivergenceL2(node,elem,soln)
    value = sqrt(integratedivergence(node,elem,soln,true));
end

function value = getdivergenceintegral(node,elem,soln)
    value = integratedivergence(node,elem,soln,false);
end

function value = integratedivergence(node,elem,soln,squareValue)
    [elem2dof,~] = dofP2(elem);
    [Dlambda,area] = gradbasis(node,elem);
    [lambda,w] = quadpts(6);
    NT = size(elem,1);
    value = 0;
    for q = 1:size(lambda,1)
        Dphi = p2gradientlocal(lambda(q,:),Dlambda);
        divq = sum(Dphi(:,1,:).*...
            reshape(soln.ux(elem2dof),NT,1,6),3) ...
            + sum(Dphi(:,2,:).*...
            reshape(soln.uz(elem2dof),NT,1,6),3);
        if squareValue
            divq = divq.^2;
        end
        value = value+w(q)*sum(area.*divq);
    end
end

function Dphi = p2gradientlocal(lambda,Dlambda)
    NT = size(Dlambda,1);
    Dphi = zeros(NT,2,6);
    Dphi(:,:,1) = (4*lambda(1)-1).*Dlambda(:,:,1);
    Dphi(:,:,2) = (4*lambda(2)-1).*Dlambda(:,:,2);
    Dphi(:,:,3) = (4*lambda(3)-1).*Dlambda(:,:,3);
    Dphi(:,:,4) = 4*(lambda(2)*Dlambda(:,:,3)+...
        lambda(3)*Dlambda(:,:,2));
    Dphi(:,:,5) = 4*(lambda(3)*Dlambda(:,:,1)+...
        lambda(1)*Dlambda(:,:,3));
    Dphi(:,:,6) = 4*(lambda(1)*Dlambda(:,:,2)+...
        lambda(2)*Dlambda(:,:,1));
end

function pressureMass = getpressuremass(node,elem)
    [~,area] = gradbasis(node,elem);
    pressureMass = accumarray(double(elem(:)),...
        repmat(area/3,3,1),[size(node,1),1]);
end
