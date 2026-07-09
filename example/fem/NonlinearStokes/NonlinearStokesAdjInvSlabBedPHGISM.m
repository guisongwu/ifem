%% NONLINEARSTOKESADJINVSLABBEDPHGISM PHGISM-parameter slab-bed beta inversion.
%
% This is a physical-unit counterpart of NonlinearStokesAdjInvSlabBed.m.
% The original script uses nondimensional data.  This version follows the
% constants in ~/software/phgism/ice-sheet/src/parameters.h, with the
% PHGISM length scaling from ins.h:
%
%   LEN_SCALING = 1e3 m, RHO_ICE = 917 kg/m^3, GRAVITY = 9.81 m/s^2,
%   POWER_N = 3, SEC_PER_YEAR = 31556926 s/year.
%
% Coordinates are in m and velocities are in m/year, matching the physical
% iFEM ISMIP-HOM examples in this directory.  The Arrhenius constants from
% parameters.h are converted to A in year^-1 Pa^-3.  PHGISM's ISMIP-HOM
% testcase fixes A=1e-16 instead; use that value below if benchmark
% matching is desired.
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

%% PHGISM constants and geometry
param = phgismparameters();

% The nondimensional slab L=4, H=1, h=0.1 is interpreted in PHGISM's
% kilometer coordinate scale and converted to meters for this solver.
L = 4*param.LEN_SCALING;
H = 1*param.LEN_SCALING;
h = 0.1*param.LEN_SCALING;
slope = 0.1;

fprintf(['PHGISM-parameter slab-bed beta inversion case: ',...
    'L = %.04e m, H = %.04e m, h = %.04e m\n'],...
    L,H,h);

[node,elem] = squaremesh([0,L,0,H],h);
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = H-slope*uNode(:,1);
tolGeometry = 100*eps(max(1,max(abs(node(:)))));
% Use all P2 velocity dofs on the top boundary as observations.  The right
% endpoint is skipped because it is identified with x=0 by periodicity.
topDof = find(abs(uNode(:,2)-surfaceLevel)<tolGeometry ...
            & uNode(:,1)<L-tolGeometry);
[~,order] = sort(uNode(topDof,1));
topDof = topDof(order);
xObs = uNode(topDof,1);
topWeight = boundaryweights(xObs,L,slope);

%% Nonlinear Stokes model
pde = struct;
pde.A = param.A;
pde.n = param.POWER_N;
pde.m = 1;
% pde.m = 1/3;
pde.rho = param.RHO_ICE;
pde.gravity = [0,-param.GRAVITY];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = param.MIN_EFFECTIVE_STRAIN*param.SEC_PER_YEAR;
option.maxIt = 200;
option.tol = 1e-6;
option.residual_tol = 1e-6;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
% The adjoint and incremental equations use the consistent nonlinear
% Jacobian assembled after the Picard forward solve.
option.assemble_tangent = true;

%% Periodic P1 beta parameterization
Nm = round(L/h);
assert(abs(Nm*h-L) <= 100*eps(max(1,L)),...
    'The parameter grid requires L/h to be an integer.');
xBeta = (0:Nm-1)'*h;
xi = mod(xBeta,L)/L;
if 0
    betaTrue = 1000*ones(size(xi));
    betaTrueName = 'constant';
elseif 0
    betaTrue = 1000*(0.85+0.30*xi);
    betaTrueName = 'linear';
elseif 0
    betaTrue = 1000*(0.9+0.4*(2*xi-1).^2);
    betaTrueName = 'quadratic';
elseif 1
    % PHGISM ISMIP-HOM-D uses 1000+1000*sin(2*pi*x/L).  The log-beta
    % parameterization used here requires beta > 0, so keep the same scale
    % and phase but use a 90% amplitude.
    betaTrue = 1000*(1+0.9*sin(2*pi*xi));
    betaTrueName = 'positive PHGISM ISMIP-HOM-D sinusoidal';
elseif 0
    betaTrue = 1000*(1+0.20*cos(2*pi*xi)+0.10*sin(4*pi*xi));
    betaTrueName = 'mixed trigonometric';
end

if 0
    betaPerturbation = 0.15*ones(size(xi));
    perturbationName = 'constant';
elseif 0
    betaPerturbation = 0.05+0.15*xi;
    perturbationName = 'linear';
elseif 0
    betaPerturbation = 0.20*((2*xi-1).^2-1/3)+0.05;
    perturbationName = 'quadratic';
elseif 1
    betaPerturbation = 0.20*sin(2*pi*xi)+0.05;
    perturbationName = 'trigonometric';
elseif 0
    betaPerturbation = 0.15*sin(2*pi*xi)+0.08*cos(4*pi*xi)+0.05;
    perturbationName = 'mixed trigonometric';
end

betaInitial = betaTrue.*(1+betaPerturbation);
qTrue = log(betaTrue);
q = log(betaInitial);
betaPlotX = linspace(0,L,401)';

fprintf(['  beta true: %s, perturbation: %s, ',...
    'beta parameters: %d, top observations: %d\n'],...
    betaTrueName,perturbationName,Nm,numel(topDof));

%% Synthetic surface observation
% This is an inverse-crime style check: generate exact data from qTrue on
% the same mesh, then try to recover q from only top-boundary velocities.
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

%% Inverse options
maxInverseIt = 20;
lambda = 1e-7;
pcgTolerance = 1e-8;
pcgMaxIt = 50;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;
useLineSearch = defaultlinesearch();

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.parameterErrorLinf = NaN(maxInverseIt,1);
history.parameterErrorRelativeLinf = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);
history.picardSteps = NaN(maxInverseIt,1);
derivativeCheck = struct('stateError',NaN,'gradientError',NaN,...
    'gaussNewtonError',NaN,'finiteDifference',NaN,...
    'adjointDirection',NaN,'forwardSolves',0);

uWarm = [];
optimizationForwardSolves = 0;
printiterationheader();
for k = 1:maxInverseIt
    iterationTimer = tic;
    % Main nonlinear state solve for the current parameter.  uWarm carries
    % the last accepted state and usually reduces the Picard iteration count.
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    optimizationForwardSolves = optimizationForwardSolves+1;
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    objective = dataObjective;

    % G is R_q, the derivative of the nonlinear residual with respect to
    % q=log(beta).  The chain rule delta beta = beta * delta q is applied
    % inside assembleparameterderivative.
    G = assembleparameterderivative(eqn,q,xBeta,L,Nm);

    % J_U is nonzero only at observed top velocity dofs.  The adjoint solve
    % uses the transposed consistent tangent: R_U' * adjoint = -J_U.
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

    if k == 1
        derivativeCheck = verifyderivatives(...
            q,u,eqn,G,gradient,dataObs,topWeight,dataNormSquared,...
            pde,option,node,elem,bdFlag,xBeta,L,topDof);
        optimizationForwardSolves = optimizationForwardSolves+...
            derivativeCheck.forwardSolves;
    end

    if objective < 1e-15
        iterationTime = toc(iterationTimer);
        printiterationrow(k,objective,history.parameterErrorRelativeLinf(k),...
            norm(gradient),forwardInfo.itStep,NaN,NaN,0,...
            iterationTime,'obj');
        history = trimhistory(history,k);
        break
    end

    if norm(gradient) <= gradientTolerance
        iterationTime = toc(iterationTimer);
        printiterationrow(k,objective,history.parameterErrorRelativeLinf(k),...
            norm(gradient),forwardInfo.itStep,NaN,NaN,0,...
            iterationTime,'grad');
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    % PCG only needs Hessian-vector products.  gaussnewtonproduct applies
    % the matrix-free GN Hessian plus the LM damping lambda*I.
    [step,flag,relativeResidual,pcgIt] = pcg(...
        hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:NonlinearAdjointPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        iterationTime = toc(iterationTimer);
        printiterationrow(k,objective,history.parameterErrorRelativeLinf(k),...
            norm(gradient),forwardInfo.itStep,pcgIt,relativeResidual,0,...
            iterationTime,'step');
        history = trimhistory(history,k);
        break
    end

    lineSearchCount = 0;
    if useLineSearch
        accepted = false;
        stepLength = 1;
        for lineSearchIt = 1:10
            lineSearchCount = lineSearchCount+1;
            qTrial = q+stepLength*step;
            % Every trial step requires a nonlinear forward solve.  These solves
            % are counted separately in optimizationForwardSolves.
            [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
                node,elem,bdFlag,xBeta,L);
            optimizationForwardSolves = optimizationForwardSolves+1;
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
    else
        qTrial = q+step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,L);
        optimizationForwardSolves = optimizationForwardSolves+1;
        if trialInfo.converged
            q = qTrial;
            uWarm = uTrial;
            lambda = max(lambda/3,1e-12);
        else
            lambda = 10*lambda;
            warning('iFEM:NonlinearAdjointNoStep',...
                'Undamped Gauss-Newton trial did not converge.');
        end
    end

    iterationTime = toc(iterationTimer);
    printiterationrow(k,objective,history.parameterErrorRelativeLinf(k),...
        norm(gradient),forwardInfo.itStep,pcgIt,relativeResidual,...
        lineSearchCount,iterationTime,'');

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = exp(q);
betaErrorLinf = norm(betaRecovered-betaTrue,inf);
betaErrorRelativeLinf = betaErrorLinf/norm(betaTrue,inf);
fprintf('\nSummary\n');
fprintf('  optimization forward solves: %d\n',...
    optimizationForwardSolves);
fprintf('  final beta Linf error: %.04e, relative Linf: %.04e\n',...
    betaErrorLinf,betaErrorRelativeLinf);
fprintf(['  derivative check: state %.04e, grad %.04e, GN %.04e ',...
         '(FD %.04e, adj %.04e)\n'],...
    derivativeCheck.stateError,derivativeCheck.gradientError,...
    derivativeCheck.gaussNewtonError,derivativeCheck.finiteDifference,...
    derivativeCheck.adjointDirection);

figure(1);
set(gcf,'Visible','on');
betaTruePlot = periodicP1(betaPlotX,xBeta,betaTrue,L);
betaInitialPlot = periodicP1(betaPlotX,xBeta,betaInitial,L);
betaRecoveredPlot = periodicP1(betaPlotX,xBeta,betaRecovered,L);
plot(betaPlotX,betaTruePlot,'k-',...
    'LineWidth',1.8,'DisplayName','true');
hold on;
plot(betaPlotX,betaInitialPlot,'b--',...
    'LineWidth',1.2,'DisplayName','initial');
plot(betaPlotX,betaRecoveredPlot,'r-',...
    'LineWidth',1.4,'DisplayName','boundary recovered');
plot(xBeta,betaRecovered,'ro','MarkerSize',5,...
    'DisplayName','recovered nodes');
hold off;
grid on;
xlabel('x (m)');
ylabel('\beta (Pa yr m^{-1})');
legend('Location','best');
title(sprintf('PHGISM-parameter slab-bed beta inversion: true %s, perturbation %s',...
    betaTrueName,perturbationName));

figure(2);
set(gcf,'Visible','on');
iteration = 1:numel(history.objective);
subplot(1,2,1);
semilogy(iteration,history.objective,'o-',...
    'LineWidth',1.4,'DisplayName','objective');
grid on;
xlabel('inverse iteration');
ylabel('objective');
legend('Location','best');
title('objective history');

iteration = 1:numel(history.parameterError);
subplot(1,2,2);
semilogy(iteration,history.parameterErrorLinf,'s-',...
    'LineWidth',1.4,'DisplayName','absolute Linf');
hold on;
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
xlabel('x (m)');
ylabel('surface horizontal velocity (m/yr)');
legend('observation','recovered prediction','Location','best');

figure(4);
set(gcf,'Visible','on');
velocityXNode = uRecovered(1:N);
velocityYNode = uRecovered(Nu+(1:N));

subplot(2,2,1);
showmesh(node,elem);
title('slab-bed mesh','FontSize',14);

subplot(2,2,2);
trisurf(elem,node(:,1),node(:,2),velocityXNode,...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered u_x','FontSize',14);
view(2);

subplot(2,2,3);
trisurf(elem,node(:,1),node(:,2),velocityYNode,...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered u_z','FontSize',14);
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

function printiterationheader()
    width = [3,12,11,12,7,5,10,2,10,4];
    label = {'it','objective','betaLinfRel','|grad|','fPicard',...
        'pcgIt','pcgRel','ls','time','stop'};
    fprintf('\n');
    fprintf('%s %s %s %s %s %s %s %s %s %s\n',...
        centertext(label{1},width(1)),centertext(label{2},width(2)),...
        centertext(label{3},width(3)),centertext(label{4},width(4)),...
        centertext(label{5},width(5)),centertext(label{6},width(6)),...
        centertext(label{7},width(7)),centertext(label{8},width(8)),...
        centertext(label{9},width(9)),centertext(label{10},width(10)));
    fprintf('%s %s %s %s %s %s %s %s %s %s\n',...
        repmat('-',1,width(1)),repmat('-',1,width(2)),...
        repmat('-',1,width(3)),repmat('-',1,width(4)),...
        repmat('-',1,width(5)),repmat('-',1,width(6)),...
        repmat('-',1,width(7)),repmat('-',1,width(8)),...
        repmat('-',1,width(9)),repmat('-',1,width(10)));
end

function printiterationrow(k,objective,betaLinfRel,...
        gradientNorm,forwardPicard,pcgIt,pcgRel,...
        lineSearchCount,iterationTime,stopReason)
    width = [3,12,11,12,7,5,10,2,10,4];
    if isnan(pcgIt)
        pcgItText = '-';
    else
        pcgItText = sprintf('%d',pcgIt);
    end
    if isnan(pcgRel)
        pcgRelText = '-';
    else
        pcgRelText = sprintf('%.04e',pcgRel);
    end
    if isempty(stopReason)
        stopReason = '-';
    end
    value = {sprintf('%d',k),sprintf('%.04e',objective),...
        sprintf('%.04e',betaLinfRel),sprintf('%.04e',gradientNorm),...
        sprintf('%d',forwardPicard),pcgItText,pcgRelText,...
        sprintf('%d',lineSearchCount),sprintf('%.2f',iterationTime),...
        stopReason};
    fprintf('%s %s %s %s %s %s %s %s %s %s\n',...
        centertext(value{1},width(1)),centertext(value{2},width(2)),...
        centertext(value{3},width(3)),centertext(value{4},width(4)),...
        centertext(value{5},width(5)),centertext(value{6},width(6)),...
        centertext(value{7},width(7)),centertext(value{8},width(8)),...
        centertext(value{9},width(9)),centertext(value{10},width(10)));
end

function text = centertext(text,width)
    text = char(text);
    pad = max(width-length(text),0);
    text = [repmat(' ',1,ceil(pad/2)),text,repmat(' ',1,floor(pad/2))];
end

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,lambda)
    % Linearized state equation:
    %     R_U * deltaU = -R_q * direction.
    incrementalState = eqn.tangent\(-G*direction);

    % Apply the observation Hessian J_uu to deltaU.  Only top-boundary
    % observed velocity components contribute to the objective.
    incrementalObservation = zeros(size(eqn.tangent,1),1);
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;

    % Linearized adjoint equation gives J_obs' * J_obs * direction without
    % assembling the dense observation Jacobian.
    incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
    product = G'*incrementalAdjoint+lambda*direction;
end

function check = verifyderivatives(q,u,eqn,G,gradient,dataObs,topWeight,...
        dataNormSquared,pde,option,node,elem,bdFlag,xBeta,L,topDof)
    % A fixed smooth direction gives a deterministic regression-style check
    % of the tangent equation, adjoint gradient, and GN quadratic form.
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

    check.stateError = relativeStateError;
    check.gradientError = relativeGradientError;
    check.gaussNewtonError = relativeGaussNewtonError;
    check.finiteDifference = finiteDifference;
    check.adjointDirection = adjointDirection;
    check.forwardSolves = 2;
end

function G = assembleparameterderivative(eqn,q,xBeta,L,Nm)
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        % Since beta=exp(q), a unit perturbation in q_j produces
        % delta beta_j = beta_j.
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
    % u0 is a warm start for the nonlinear Picard iteration; omit it for
    % truth solves or the first inverse iteration.
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
        history.(fields{j}) = history.(fields{j})(1:k,:);
    end
end

function value = defaultlinesearch()
    value = false;
end

function param = phgismparameters()
    param.ARRHENIUS_T = 263.15;
    param.ARRHENIUS_A0 = 3.61e-13;
    param.ARRHENIUS_A1 = 1.73e3;
    param.ARRHENIUS_Q0 = 6.0e4;
    param.ARRHENIUS_Q1 = 13.9e4;
    param.RHO_ICE = 917;
    param.RHO_WATER = 1027;
    param.GRAVITY = 9.81;
    param.POWER_N = 3;
    param.TEMP_WATER = 273.15;
    param.GEOTHE_FLUX = 3e-2;
    param.THEM_CONDUCT = 2.1;
    param.HEAT_CAPACITY = 2009;
    param.BETA_MELT = 8.66e-4;
    param.LATENT_CAPACITY = 3.35e5;
    param.GAS_CONSTANT = 8.314;
    param.SEC_PER_YEAR = 31556926;
    param.MIN_EFFECTIVE_STRAIN = 1e-15;
    param.HEAT_SOURCE = 0;

    % Defined in PHGISM ins.h and used by HEIGHT_EPS in parameters.h.
    param.LEN_SCALING = 1e3;
    param.HEIGHT_EPS = 1/param.LEN_SCALING;

    if param.ARRHENIUS_T < 263.15
        flowRatePerSecond = param.ARRHENIUS_A0*exp(...
            -param.ARRHENIUS_Q0/(param.GAS_CONSTANT*param.ARRHENIUS_T));
    else
        flowRatePerSecond = param.ARRHENIUS_A1*exp(...
            -param.ARRHENIUS_Q1/(param.GAS_CONSTANT*param.ARRHENIUS_T));
    end
    param.A = flowRatePerSecond*param.SEC_PER_YEAR;
end
