%% FULLSTOKESADJINVSINISM PHGISM-parameter sinusoidal-bed beta inversion.
%
% This is the PHGISM-internal-unit counterpart of
% FullStokesAdjInvSin.m.  The original script uses nondimensional
% data on an ISMIP-HOM-B-like sinusoidal bed.  This version follows the constants in
% ~/software/phgism/ice-sheet/src/parameters.h and the scaling macros in
% ~/software/phgism/ice-sheet/src/ins.h:
%
%   LEN_SCALING = 1e3 m, RHO_ICE = 917 kg/m^3, GRAVITY = 9.81 m/s^2,
%   POWER_N = 3, SEC_PER_YEAR = 31556926 s/year,
%   EQU_SCALING = 1e-8, PRES_SCALING = 1e5 Pa.
%
% Coordinates are in km, velocities are in m/year, and returned pressure is
% the PHGISM pressure degree of freedom p/1e5.  Because NonlinearStokesP2P1
% differentiates with respect to the supplied coordinates and uses an
% unscaled pressure block, the coefficients passed to it are the PHGISM
% internal Stokes coefficients:
%
%   A_ifem        = A / (LEN_SCALING^(n-1) * EQU_SCALING^n),
%   eps_reg_ifem  = eps_reg_physical * LEN_SCALING,
%   f_ifem        = rho*g * LEN_SCALING^2 * EQU_SCALING,
%   beta_ifem     = beta * LEN_SCALING * EQU_SCALING.
%
% The iFEM pressure unknown is therefore
% EQU_SCALING*LEN_SCALING*PRES_SCALING times the PHGISM pressure DOF.
% The Arrhenius constants from parameters.h are converted to A in
% year^-1 Pa^-3 before the PHGISM matrix scaling is applied.
%
%     0.5 * int_{\Gamma_t} (u-u_obs)^2 ds
%         / int_{\Gamma_t} u_obs^2 ds.
%
% A periodic first-difference Tikhonov term can be added to the inverse
% objective through alpha.  The default alpha=0 keeps the data-only case.
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

% The nondimensional sinusoidal-bed case is interpreted directly in
% PHGISM's internal kilometer coordinate scale.
L = 4;
H = 1;
h = 0.1;
slope = tan(0.5*pi/180);
bedAmplitude = 0.1*H;
Nx = max(4,round(L/h));
Nz = max(3,round(H/h));

fprintf('ISM sin-bed inversion: L=%.2f km, H=%.2f km, h=%.2f km, Nx=%d, Nz=%d\n',...
    L,H,h,Nx,Nz);

[node,elem] = rectanglemesh(L,1,Nx,Nz);
% Boundary flags are set on the reference strip [0,L]x[0,1] before the
% B-type geometry map.  The physical top is z=-slope*x after mapping, but
% its reference location is always z=1.
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node = maptoexperimentb(node,L,H,bedAmplitude,slope);

[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = -slope*uNode(:,1);
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
pde.A = param.A/(param.LEN_SCALING^(param.POWER_N-1) * ...
    param.EQU_SCALING^param.POWER_N);
pde.n = param.POWER_N;
pde.m = 1;
% pde.m = 1/3;
pde.rho = param.RHO_ICE;
pde.gravity = [0,-param.GRAVITY*param.LEN_SCALING^2*param.EQU_SCALING];
pde.g_N = [];
pde.beta_scale = param.LEN_SCALING*param.EQU_SCALING;
pde.pressure_dof_scale = param.EQU_SCALING*param.LEN_SCALING*...
    param.PRES_SCALING;

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = param.MIN_EFFECTIVE_STRAIN*param.LEN_SCALING;
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

% Choose one positive sinusoidal beta truth and one sinusoidal initial
% perturbation.  betaTrue must remain positive because q = log(beta).

% PHGISM ISMIP-HOM-D shape with 90% amplitude, range [100, 1900].
betaTrueCase = 'homd90';
betaTrue = 1000*(1+0.9*sin(2*pi*xi));
betaTrueName = 'sin 1 wave amp 0.90 phase 0';


% betaPerturbationCase = 'sin20_offset05';
% betaPerturbationCase = 'sin10_offset05';
% betaPerturbationCase = 'shifted15_offset05';
betaPerturbationCase = 'doublewave10_offset05';
% betaPerturbationCase = 'opposite20_offset05';


switch lower(betaPerturbationCase)
    case 'sin20_offset05'
        betaPerturbation = 0.05+0.20*sin(2*pi*xi);
        perturbationName = 'sin 1 wave amp 0.20 offset 0.05';
    case 'sin10_offset05'
        betaPerturbation = 0.05+0.10*sin(2*pi*xi);
        perturbationName = 'sin 1 wave amp 0.10 offset 0.05';
    case 'shifted15_offset05'
        betaPerturbation = 0.05+0.15*sin(2*pi*xi+pi/3);
        perturbationName = 'sin 1 wave amp 0.15 phase pi/3 offset 0.05';
    case 'doublewave10_offset05'
        betaPerturbation = 0.05+0.10*sin(4*pi*xi);
        perturbationName = 'sin 2 waves amp 0.10 offset 0.05';
    case 'opposite20_offset05'
        betaPerturbation = 0.05-0.20*sin(2*pi*xi);
        perturbationName = 'sin 1 wave amp -0.20 offset 0.05';
    otherwise
        error('Unknown betaPerturbationCase: %s',betaPerturbationCase);
end

betaInitial = betaTrue.*(1+betaPerturbation);
qTrue = log(betaTrue);
q = log(betaInitial);
betaPlotX = linspace(0,L,401)';

fprintf('  beta true: %s, beta init: %s\n',...
    betaTrueName,perturbationName);
fprintf('  beta parameters: %d, top obs: %d\n',...
    Nm,numel(topDof));

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
alpha = 0;
lambda = 1e-7;
pcgTolerance = 1e-8;
pcgMaxIt = 50;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;
useLineSearch = defaultlinesearch();
D = spdiags([-ones(Nm,1),ones(Nm,1)],[0,1],Nm,Nm);
D(Nm,1) = 1;

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
    regularization = D*q;
    regularizationObjective = 0.5*alpha*...
        (regularization'*regularization);
    objective = dataObjective+regularizationObjective;

    % G is R_q, the derivative of the nonlinear residual with respect to
    % q=log(beta).  The chain rule delta beta = beta * delta q is applied
    % inside assembleparameterderivative.
    G = assembleparameterderivative(eqn,q,xBeta,L,Nm,pde.beta_scale);

    % J_U is nonzero only at observed top velocity dofs.  The adjoint solve
    % uses the transposed consistent tangent: R_U' * adjoint = -J_U.
    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint+alpha*(D'*(D*q));

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
            pde,option,node,elem,bdFlag,xBeta,L,topDof,D,alpha);
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
        topDof,topWeight,dataNormSquared,D,alpha,lambda);
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
                trialRegularization = D*qTrial;
                trialObjective = 0.5*(topWeight'*(trialResidual.^2)) ...
                    /dataNormSquared+0.5*alpha*...
                    (trialRegularization'*trialRegularization);
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
xlabel('x (km)');
ylabel('\beta (Pa yr m^{-1})');
legend('Location','best');
title(sprintf('PHGISM-parameter sin-bed beta inversion: true %s, perturbation %s',...
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
xlabel('x (km)');
ylabel('surface horizontal velocity (m/yr)');
legend('observation','recovered prediction','Location','best');

figure(4);
set(gcf,'Visible','on');
velocityXNode = uRecovered(1:N);
velocityYNode = uRecovered(Nu+(1:N));

subplot(2,2,1);
showmesh(node,elem);
title('sinusoidal-bed mesh','FontSize',14);

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
title('recovered pressure DOF p/1e5','FontSize',14);
drawnow;
exportepsfigures(mfilename);


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
        topWeight,dataNormSquared,D,alpha,lambda)
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
    product = G'*incrementalAdjoint+alpha*(D'*(D*direction))+...
        lambda*direction;
end

function check = verifyderivatives(q,u,eqn,G,gradient,dataObs,topWeight,...
        dataNormSquared,pde,option,node,elem,bdFlag,xBeta,L,topDof,...
        D,alpha)
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
        /dataNormSquared+0.5*alpha*...
        ((D*(q+epsilon*direction))'*(D*(q+epsilon*direction)));
    minusObjective = 0.5*(topWeight'*(minusResidual.^2)) ...
        /dataNormSquared+0.5*alpha*...
        ((D*(q-epsilon*direction))'*(D*(q-epsilon*direction)));
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
    tangentRegularization = sqrt(alpha)*(D*direction);
    finiteDifferenceObservation = ...
        sqrt(topWeight/dataNormSquared).*...
        (uPlus(topDof)-uMinus(topDof))/(2*epsilon);
    finiteDifferenceRegularization = ...
        sqrt(alpha)*(D*(q+epsilon*direction)...
        -D*(q-epsilon*direction))/(2*epsilon);
    relativeGaussNewtonError = abs(...
        tangentObservation'*tangentObservation-...
        finiteDifferenceObservation'*finiteDifferenceObservation+...
        tangentRegularization'*tangentRegularization-...
        finiteDifferenceRegularization'*finiteDifferenceRegularization)/...
        max([eps,tangentObservation'*tangentObservation,...
             finiteDifferenceObservation'*finiteDifferenceObservation,...
             tangentRegularization'*tangentRegularization,...
             finiteDifferenceRegularization'*...
             finiteDifferenceRegularization]);

    check.stateError = relativeStateError;
    check.gradientError = relativeGradientError;
    check.gaussNewtonError = relativeGaussNewtonError;
    check.finiteDifference = finiteDifference;
    check.adjointDirection = adjointDirection;
    check.forwardSolves = 2;
end

function G = assembleparameterderivative(eqn,q,xBeta,L,Nm,betaScale)
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        % Since beta=exp(q), a unit perturbation in q_j produces
        % delta beta_j = beta_j.
        deltaBeta = betaScale*beta.*direction;
        directionFunction = @(pt) periodicP1(...
            pt(:,1),xBeta,deltaBeta,L);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L)
    beta = exp(q(:));
    pde.beta = @(pt) pde.beta_scale*periodicP1(pt(:,1),xBeta,beta,L);
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
    p = soln.p/pde.pressure_dof_scale;
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

function node = maptoexperimentb(node,L,H,bedAmplitude,slope)
    x = node(:,1);
    r = node(:,2);
    surface = -slope*x;
    bed = surface-H+bedAmplitude*sin(2*pi*x/L);
    node(:,2) = bed+r.*(surface-bed);
end

function [node,elem] = rectanglemesh(L,H,Nx,Nz)
    x = linspace(0,L,Nx+1);
    z = linspace(0,H,Nz+1);
    [X,Z] = meshgrid(x,z);
    node = [X(:),Z(:)];

    cellId = reshape(1:(Nx+1)*(Nz+1),Nz+1,Nx+1);
    elem = zeros(2*Nx*Nz,3);
    cursor = 0;
    for ix = 1:Nx
        for iz = 1:Nz
            v1 = cellId(iz,ix);
            v2 = cellId(iz,ix+1);
            v3 = cellId(iz+1,ix);
            v4 = cellId(iz+1,ix+1);
            cursor = cursor+1;
            elem(cursor,:) = [v1,v2,v4];
            cursor = cursor+1;
            elem(cursor,:) = [v1,v4,v3];
        end
    end
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
    param.EQU_SCALING = 1e-8;
    param.LEN_SCALING = 1e3;
    param.PRES_SCALING = 1e5;
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

function exportepsfigures(scriptName)
    outputDir = fullfile(fileparts(mfilename('fullpath')),'output',scriptName);
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end

    figs = findall(0,'Type','figure');
    if isempty(figs)
        return;
    end
    figNumbers = arrayfun(@(fig) fig.Number,figs);
    [~,order] = sort(figNumbers);
    figs = figs(order);

    for i = 1:numel(figs)
        fig = figs(i);
        if isgraphics(fig,'figure')
            set(fig,'Renderer','painters');
            filename = fullfile(outputDir,sprintf('%s_figure_%02d.eps',...
                scriptName,fig.Number));
            print(fig,filename,'-depsc','-vector');
        end
    end
    fprintf('  exported EPS figures to %s\n',outputDir);
end
