%% NSDERIVATIVECOMPARISON Compare FD and adjoint derivatives.
%
% This script uses the same setup as NonlinearStokesAdjInvSlabBed.m and
% compares
%
%   1. the objective gradient from centered finite differences and the
%      nonlinear adjoint equation;
%   2. the Gauss-Newton Hessian from an FD observation Jacobian and the
%      matrix-free tangent/adjoint construction.
%
% The Hessian compared here is the Gauss-Newton approximation used by the
% inversion script, not the full second derivative of the reduced objective.
% All derivatives are with respect to q=log(beta), not beta itself.

close all;
clear variables;

%% Geometry and mesh
L = 1;
H = 0.5;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

[~,edge] = dofP2(elem);
N = size(node,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = H-slope*uNode(:,1);
tolGeometry = 100*eps(max(1,max(abs(node(:)))));
topDof = find(abs(uNode(:,2)-surfaceLevel)<tolGeometry ...
            & uNode(:,1)<L-tolGeometry);
[~,order] = sort(uNode(topDof,1));
topDof = topDof(order);
xObs = uNode(topDof,1);
topWeight = boundaryweights(xObs,L,slope);

%% Nonlinear Stokes model
pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1/3;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.assemble_tangent = true;

%% Periodic P1 parameterization
Nm = round(L/h);
assert(abs(Nm*h-L) <= 100*eps(max(1,L)),...
    'The parameter grid requires L/h to be an integer.');
xBeta = (0:Nm-1)'*h;
betaTrue = 1+0.1*cos(2*pi*xBeta/L);
betaInitial = betaTrue+0.1*(sin(2*pi*xBeta/L)+0.25);
qTrue = log(betaTrue);
q = log(betaInitial);

%% Synthetic surface observation
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

%% Reference state and adjoint derivatives at q
[u,eqn,forwardInfo] = solveforward(q,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(forwardInfo.converged,'The reference solve did not converge.');

residual = u(topDof)-dataObs;
objective = objectivefromstate(u,dataObs,topDof,topWeight,...
    dataNormSquared);

G = assembleparameterderivative(eqn,q,xBeta,L,Nm);
observationGradient = zeros(size(eqn.tangent,1),1);
observationGradient(topDof) = topWeight.*residual/dataNormSquared;
adjoint = eqn.tangent'\(-observationGradient);
adjointGradient = G'*adjoint;
adjointHessianGN = assemblegaussnewtonhessian(eqn,G,topDof,...
    topWeight,dataNormSquared,Nm);

%% Finite-difference derivatives
gradientStep = 1e-4;
jacobianStep = 1e-4;
fdGradient = finiteDifferenceGradient(q,u,pde,option,node,elem,bdFlag,...
    xBeta,L,dataObs,topDof,topWeight,dataNormSquared,gradientStep);
fdJacobian = finiteDifferenceObservationJacobian(q,u,pde,option,node,...
    elem,bdFlag,xBeta,L,topDof,topWeight,dataNormSquared,jacobianStep);
fdHessianGN = fdJacobian'*fdJacobian;

%% Report
gradientDiff = fdGradient-adjointGradient;
hessianDiff = fdHessianGN-adjointHessianGN;

fprintf('\nNonlinear Stokes derivative comparison\n');
fprintf('  objective: %.04e\n',objective);
fprintf('  reference forward Picard steps: %d\n',forwardInfo.itStep);
fprintf('  gradient FD step: %.04e\n',gradientStep);
fprintf('  Jacobian FD step: %.04e\n',jacobianStep);
fprintf('\nGradient comparison\n');
fprintf('  ||g_FD||                 %.04e\n',norm(fdGradient));
fprintf('  ||g_adjoint||            %.04e\n',norm(adjointGradient));
fprintf('  ||g_FD-g_adjoint||       %.04e\n',norm(gradientDiff));
fprintf('  relative difference      %.04e\n',...
    norm(gradientDiff)/max([eps,norm(fdGradient),norm(adjointGradient)]));
fprintf('  max absolute difference  %.04e\n',norm(gradientDiff,inf));

fprintf('\nGauss-Newton Hessian comparison\n');
fprintf('  ||H_FD||_F               %.04e\n',norm(fdHessianGN,'fro'));
fprintf('  ||H_adjoint||_F          %.04e\n',norm(adjointHessianGN,'fro'));
fprintf('  ||H_FD-H_adjoint||_F     %.04e\n',norm(hessianDiff,'fro'));
fprintf('  relative difference      %.04e\n',...
    norm(hessianDiff,'fro')/max([eps,norm(fdHessianGN,'fro'),...
    norm(adjointHessianGN,'fro')]));
fprintf('  max absolute difference  %.04e\n',norm(hessianDiff(:),inf));
fprintf('  symmetry error, H_adjoint %.04e\n',...
    norm(adjointHessianGN-adjointHessianGN','fro')/...
    max(eps,norm(adjointHessianGN,'fro')));

fprintf('\nFirst few entries\n');
fprintf('  gradient columns: [FD, adjoint, difference]\n');
disp([fdGradient,adjointGradient,gradientDiff]);
fprintf('  diagonal Hessian columns: [FD, adjoint, difference]\n');
disp([diag(fdHessianGN),diag(adjointHessianGN),diag(hessianDiff)]);

%% Visualization
fprintf('\nWhat is plotted\n');
fprintf(['  Figure 1: gradient dJ/dq, analogous to dXidm in ',...
    '../Stokes/stokes_inversion.m, but for q=log(beta).\n']);
fprintf(['  Figure 2: Gauss-Newton Hessian H=J_obs''*J_obs, ',...
    'where J_obs maps q to weighted top velocity observations.\n']);
fprintf(['  Figure 3: Hessian-vector product for a smooth test direction, ',...
    'showing H*d from both constructions.\n']);

figure(1);
set(gcf,'Visible','on');
plotsensitivitycurve(xBeta,L,fdGradient,'ko-',...
    'FD gradient');
hold on;
plotsensitivitycurve(xBeta,L,adjointGradient,'r.-',...
    'adjoint gradient');
plotsensitivitycurve(xBeta,L,gradientDiff,'b.--',...
    'FD - adjoint');
hold off;
grid on;
xlabel('x');
ylabel('d objective / d q');
legend('Location','best');
title('Gradient comparison with respect to q=log(\beta)');

figure(2);
set(gcf,'Visible','on');
hessianLimit = max(abs([fdHessianGN(:);adjointHessianGN(:)]));
if hessianLimit == 0
    hessianLimit = 1;
end
diffLimit = max(abs(hessianDiff(:)));
if diffLimit == 0
    diffLimit = 1;
end
subplot(1,3,1);
imagesc(fdHessianGN);
axis image;
colorbar;
clim([-hessianLimit,hessianLimit]);
title('FD GN Hessian');
xlabel('q column');
ylabel('q row');
subplot(1,3,2);
imagesc(adjointHessianGN);
axis image;
colorbar;
clim([-hessianLimit,hessianLimit]);
title('adjoint GN Hessian');
xlabel('q column');
ylabel('q row');
subplot(1,3,3);
imagesc(hessianDiff);
axis image;
colorbar;
clim([-diffLimit,diffLimit]);
title('FD - adjoint');
xlabel('q column');
ylabel('q row');
colormap(parula);

figure(3);
set(gcf,'Visible','on');
testDirection = sin((1:Nm)');
testDirection = testDirection/norm(testDirection);
fdHessianAction = fdHessianGN*testDirection;
adjointHessianAction = adjointHessianGN*testDirection;
plotsensitivitycurve(xBeta,L,fdHessianAction,'ko-',...
    'FD H*d');
hold on;
plotsensitivitycurve(xBeta,L,adjointHessianAction,'r.-',...
    'adjoint H*d');
plotsensitivitycurve(xBeta,L,fdHessianAction-adjointHessianAction,...
    'b.--','difference');
hold off;
grid on;
xlabel('x');
ylabel('Hessian-vector product');
legend('Location','best');
title('Gauss-Newton Hessian action comparison');
drawnow;

function Hgn = assemblegaussnewtonhessian(eqn,G,topDof,topWeight,...
        dataNormSquared,Nm)
    Hgn = zeros(Nm,Nm);
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        Hgn(:,j) = gaussnewtonproduct(direction,eqn,G,topDof,...
            topWeight,dataNormSquared,0);
    end
    Hgn = 0.5*(Hgn+Hgn');
end

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,lambda)
    incrementalState = eqn.tangent\(-G*direction);
    incrementalObservation = zeros(size(eqn.tangent,1),1);
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;
    incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
    product = G'*incrementalAdjoint+lambda*direction;
end

function gradient = finiteDifferenceGradient(q,uWarm,pde,option,node,elem,...
        bdFlag,xBeta,L,dataObs,topDof,topWeight,dataNormSquared,step)
    Nm = numel(q);
    gradient = zeros(Nm,1);
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        plusObjective = objectivefromparameter(q+step*direction,uWarm,...
            pde,option,node,elem,bdFlag,xBeta,L,dataObs,topDof,...
            topWeight,dataNormSquared);
        minusObjective = objectivefromparameter(q-step*direction,uWarm,...
            pde,option,node,elem,bdFlag,xBeta,L,dataObs,topDof,...
            topWeight,dataNormSquared);
        gradient(j) = (plusObjective-minusObjective)/(2*step);
    end
end

function J = finiteDifferenceObservationJacobian(q,u,pde,option,node,elem,...
        bdFlag,xBeta,L,topDof,topWeight,dataNormSquared,step)
    Nm = numel(q);
    J = zeros(numel(topDof),Nm);
    observationWeight = sqrt(topWeight/dataNormSquared);
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        [uPlus,~,plusInfo] = solveforward(q+step*direction,u,...
            pde,option,node,elem,bdFlag,xBeta,L);
        [uMinus,~,minusInfo] = solveforward(q-step*direction,u,...
            pde,option,node,elem,bdFlag,xBeta,L);
        assert(plusInfo.converged && minusInfo.converged,...
            'Finite-difference Jacobian solve failed for column %d.',j);
        J(:,j) = observationWeight.*...
            (uPlus(topDof)-uMinus(topDof))/(2*step);
    end
end

function objective = objectivefromparameter(q,uWarm,pde,option,node,elem,...
        bdFlag,xBeta,L,dataObs,topDof,topWeight,dataNormSquared)
    [u,~,info] = solveforward(q,uWarm,pde,option,node,elem,bdFlag,...
        xBeta,L);
    assert(info.converged,'Finite-difference objective solve failed.');
    objective = objectivefromstate(u,dataObs,topDof,topWeight,...
        dataNormSquared);
end

function objective = objectivefromstate(u,dataObs,topDof,topWeight,...
        dataNormSquared)
    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
end

function G = assembleparameterderivative(eqn,q,xBeta,L,Nm)
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        deltaBeta = beta.*direction;
        directionFunction = @(pt) periodicP1(...
            pt(:,1),xBeta,deltaBeta,L);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L)
    beta = exp(q(:));
    pde.beta = @(pt) periodicP1(pt(:,1),xBeta,beta,L);
    if isempty(u0)
        if isfield(option,'u0')
            option = rmfield(option,'u0');
        end
    else
        option.u0 = u0;
    end
    [soln,eqn,info] = NonlinearStokesP2P1(...
        node,elem,bdFlag,pde,option);
    u = soln.u;
    p = soln.p;
end

function value = periodicP1(x,xNode,nodalValue,L)
    xWrapped = mod(x,L);
    value = interp1([xNode;L],[nodalValue;nodalValue(1)],...
        xWrapped,'linear');
end

function weight = boundaryweights(xObs,L,slope)
    nObs = numel(xObs);
    assert(nObs > 0,'No top-boundary observation dofs were found.');
    weight = sqrt(1+slope^2)*(L/nObs)*ones(nObs,1);
end

function plotsensitivitycurve(xBeta,L,value,lineSpec,displayName)
    plot([xBeta;L],[value;value(1)],lineSpec,...
        'LineWidth',1.3,'MarkerSize',5,'DisplayName',displayName);
end
