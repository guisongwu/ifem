%% NONLINEARSTOKESADJOINTINVERSIONBOUNDARYOBJECTIVE Adjoint inversion with boundary-integral misfit.
%
% This is a copy of NonlinearStokesAdjointInversion.m with the data
% objective changed to a top-boundary integral
%
%     0.5 * int_{\Gamma_t} (u-u_obs)^2 ds
%         / int_{\Gamma_t} u_obs^2 ds.
%
% No regularization term is included in the inverse objective.
%
% The state equation uses nonlinear Glen viscosity and a regularized
% Weertman sliding law.  The inversion variable is q=log(beta).
% Gradients are computed by one consistent nonlinear adjoint solve.
% Gauss-Newton steps use matrix-free incremental-state/adjoint solves.

close all;
clear variables;
set(groot,'DefaultFigureVisible','on');

%% Geometry and mesh
L = 1;
H = 0.5;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==0.5','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1); %#ok<NASGU>
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

%% Inverse options
maxInverseIt = 10;
lambda = 1e-7;
pcgTolerance = 1e-8;
pcgMaxIt = 50;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.parameterErrorLinf = NaN(maxInverseIt,1);
history.parameterErrorRelativeLinf = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);
history.picardSteps = NaN(maxInverseIt,1);
derivativeCheck = struct('stateError',NaN,'gradientError',NaN,...
    'gaussNewtonError',NaN,'finiteDifference',NaN,...
    'adjointDirection',NaN);

uWarm = [];
for k = 1:maxInverseIt
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    objective = dataObjective;

    G = assembleparameterderivative(eqn,q,xBeta,L,Nm);
    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

    history.objective(k) = objective;
    history.dataResidual(k) = sqrt(...
        (topWeight'*(residual.^2))/dataNormSquared);
    betaCurrent = exp(q);
    history.parameterError(k) = norm(betaCurrent-betaTrue)/norm(betaTrue);
    history.parameterErrorLinf(k) = norm(betaCurrent-betaTrue,inf);
    history.parameterErrorRelativeLinf(k) = ...
        history.parameterErrorLinf(k)/norm(betaTrue,inf);
    history.gradientNorm(k) = norm(gradient);
    history.picardSteps(k) = forwardInfo.itStep;

    fprintf(['boundary inverse %2d: objective %.6e, ',...
             'beta L2 rel %.6e, beta Linf abs %.6e, ',...
             'beta Linf rel %.6e, |gradient| %.6e, Picard %d\n'],...
        k,objective,history.parameterError(k),...
        history.parameterErrorLinf(k),...
        history.parameterErrorRelativeLinf(k),norm(gradient),...
        forwardInfo.itStep);

    if k == 1
        derivativeCheck = verifyderivatives(...
            q,u,eqn,G,gradient,dataObs,topWeight,dataNormSquared,...
            pde,option,node,elem,bdFlag,xBeta,L,topDof);
    end

    if norm(gradient) <= gradientTolerance
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    [step,flag,relativeResidual,pcgIt] = pcg(...
        hessian,-gradient,pcgTolerance,pcgMaxIt);
    fprintf('  PCG flag %d, iteration %d, relative residual %.3e\n',...
        flag,pcgIt,relativeResidual);
    if flag ~= 0
        warning('iFEM:NonlinearAdjointPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        history = trimhistory(history,k);
        break
    end

    accepted = false;
    stepLength = 1;
    for lineSearchIt = 1:10
        qTrial = q+stepLength*step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,L);
        if trialInfo.converged
            trialResidual = uTrial(topDof)-dataObs;
            trialObjective = 0.5*(topWeight'*(trialResidual.^2)) ...
                /dataNormSquared;
            if trialObjective < objective
                q = qTrial;
                uWarm = uTrial;
                lambda = max(lambda/3,1e-12);
                accepted = true;
                break
            end
        end
        stepLength = stepLength/2;
    end

    if ~accepted
        lambda = 10*lambda;
        warning('iFEM:NonlinearAdjointNoStep',...
            'No decreasing step was found; increasing LM damping.');
    end

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = exp(q);
betaErrorLinf = norm(betaRecovered-betaTrue,inf);
betaErrorRelativeLinf = betaErrorLinf/norm(betaTrue,inf);
fprintf('final beta Linf error %.6e, relative Linf %.6e\n',...
    betaErrorLinf,betaErrorRelativeLinf);

figure(1);
set(gcf,'Visible','on');
plotFine = linspace(0,L,401)';
betaTruePlot = periodicP1(plotFine,xBeta,betaTrue,L);
betaInitialPlot = periodicP1(plotFine,xBeta,betaInitial,L);
betaRecoveredPlot = periodicP1(plotFine,xBeta,betaRecovered,L);
plot(plotFine,betaTruePlot,'k-',...
    'LineWidth',1.8,'DisplayName','true');
hold on;
plot(plotFine,betaInitialPlot,'b--',...
    'LineWidth',1.2,'DisplayName','initial');
plot(plotFine,betaRecoveredPlot,'r-',...
    'LineWidth',1.4,'DisplayName','boundary recovered');
plot(xBeta,betaRecovered,'ro','MarkerSize',5,...
    'DisplayName','recovered nodes');
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title('Nonlinear Stokes adjoint inversion, boundary objective');

figure(2);
set(gcf,'Visible','on');
iteration = 1:numel(history.parameterError);
semilogy(iteration,history.parameterError,'o-',...
    'LineWidth',1.4,'DisplayName','relative L2');
hold on;
semilogy(iteration,history.parameterErrorLinf,'s-',...
    'LineWidth',1.4,'DisplayName','absolute Linf');
semilogy(iteration,history.parameterErrorRelativeLinf,'^-',...
    'LineWidth',1.4,'DisplayName','relative Linf');
hold off;
grid on;
xlabel('inverse iteration');
ylabel('\beta error');
legend('Location','best');
title('\beta error history');

figure(3);
set(gcf,'Visible','on');
[uRecovered,~,~,pRecovered] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,L);
plot(xObs,dataObs,'ko',xObs,uRecovered(topDof),'r-',...
    'LineWidth',1.4);
grid on;
xlabel('x');
ylabel('surface horizontal velocity');
legend('observation','recovered prediction','Location','best');

figure(4);
set(gcf,'Visible','on');
velocityXNode = uRecovered(1:N);
velocityYNode = uRecovered(Nu+(1:N));
velocityMagnitudeNode = sqrt(velocityXNode.^2+velocityYNode.^2);

subplot(2,2,1);
trisurf(elem,node(:,1),node(:,2),velocityXNode,...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered u_x','FontSize',14);
view(2);

subplot(2,2,2);
trisurf(elem,node(:,1),node(:,2),velocityYNode,...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered u_y','FontSize',14);
view(2);

subplot(2,2,3);
trisurf(elem,node(:,1),node(:,2),velocityMagnitudeNode,...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered |u|','FontSize',14);
view(2);

subplot(2,2,4);
trisurf(elem,node(:,1),node(:,2),pRecovered,...
    'FaceColor','interp','EdgeColor','interp');
view(2);
axis equal;
axis tight;
colorbar;
title('recovered pressure','FontSize',14);
drawnow;

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,lambda)
    incrementalState = eqn.tangent\(-G*direction);
    incrementalObservation = zeros(size(eqn.tangent,1),1);
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;
    incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
    product = G'*incrementalAdjoint+lambda*direction;
end

function check = verifyderivatives(q,u,eqn,G,gradient,dataObs,topWeight,...
        dataNormSquared,pde,option,node,elem,bdFlag,xBeta,L,topDof)
    direction = sin((1:numel(q))');
    direction = direction/norm(direction);
    epsilon = 1e-3;

    [uPlus,~,plusInfo] = solveforward(q+epsilon*direction,u,...
        pde,option,node,elem,bdFlag,xBeta,L);
    [uMinus,~,minusInfo] = solveforward(q-epsilon*direction,u,...
        pde,option,node,elem,bdFlag,xBeta,L);
    assert(plusInfo.converged && minusInfo.converged,...
        'Gradient-check forward solve failed.');
    plusResidual = uPlus(topDof)-dataObs;
    minusResidual = uMinus(topDof)-dataObs;
    plusObjective = 0.5*(topWeight'*(plusResidual.^2)) ...
        /dataNormSquared;
    minusObjective = 0.5*(topWeight'*(minusResidual.^2)) ...
        /dataNormSquared;
    finiteDifference = (plusObjective-minusObjective)/(2*epsilon);
    adjointDirection = gradient'*direction;
    relativeGradientError = abs(finiteDifference-adjointDirection)/...
        max([eps,abs(finiteDifference),abs(adjointDirection)]);

    incrementalState = eqn.tangent\(-G*direction);
    finiteDifferenceState = (uPlus-uMinus)/(2*epsilon);
    relativeStateError = norm(...
        finiteDifferenceState-incrementalState(1:numel(u)))/...
        max(eps,norm(finiteDifferenceState));

    tangentObservation = ...
        sqrt(topWeight/dataNormSquared).*incrementalState(topDof);
    finiteDifferenceObservation = ...
        sqrt(topWeight/dataNormSquared).*...
        (uPlus(topDof)-uMinus(topDof))/(2*epsilon);
    relativeGaussNewtonError = abs(...
        tangentObservation'*tangentObservation-...
        finiteDifferenceObservation'*finiteDifferenceObservation)/...
        max([eps,tangentObservation'*tangentObservation,...
             finiteDifferenceObservation'*finiteDifferenceObservation]);

    fprintf(['  derivative check: state %.3e, gradient %.3e, ',...
             'Gauss-Newton %.3e ',...
             '(FD %.6e, adjoint %.6e)\n'],...
        relativeStateError,relativeGradientError,...
        relativeGaussNewtonError,...
        finiteDifference,adjointDirection);
    check.stateError = relativeStateError;
    check.gradientError = relativeGradientError;
    check.gaussNewtonError = relativeGaussNewtonError;
    check.finiteDifference = finiteDifference;
    check.adjointDirection = adjointDirection;
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

function history = trimhistory(history,k)
    fields = fieldnames(history);
    for j = 1:numel(fields)
        history.(fields{j}) = history.(fields{j})(1:k);
    end
end
