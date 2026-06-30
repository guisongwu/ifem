%% FOADJINVSINBED Adjoint beta inversion for the 2D FO sinusoidal bed.
%
% First-order counterpart of ../NonlinearStokes/NonlinearStokesAdjInvSinBed.m.
% Synthetic top-surface horizontal velocity observations are generated on
% the same mapped geometry and used to recover beta on the basal boundary.

if exist('foInvConfig','var')
    close all;
else
    close all;
    clear variables;
    foInvConfig = struct;
end
set(groot,'DefaultFigureVisible',getconfig(foInvConfig,'figureVisible','on'));

%% Geometry and model
% Fixed-thickness cases.  Toggle exactly one branch to compare how
% horizontal length affects beta identifiability from top-surface velocity
% data.
if 0
    L = 1;
    domainCase = 'L1';
elseif 0
    L = 2;
    domainCase = 'L2';
elseif 1
    L = 4;
    domainCase = 'L4';
end
H = 1;
h = 0.1;
slope = tan(0.5*pi/180);
bedAmplitude = 0.1*H;
Nx = max(4,round(L/h));
Nz = max(3,round(H/h));

fprintf(['FO sin-bed beta inversion case %s: L = %.04e, ',...
    'H = %.04e, h = %.04e, Nx = %d, Nz = %d\n'],...
    domainCase,L,H,h,Nx,Nz);

[node,elem] = rectanglemesh(L,1,Nx,Nz);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node = maptoexperimentb(node,L,H,bedAmplitude,slope);

pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = 1;
pde.f = [];

option.L = L;
option.H = H;
option.slope = slope;
option.node = node;
option.elem = elem;
option.bdFlag = bdFlag;
option.periodic_x = [0,L];
option.A = pde.A;
option.n = pde.n;
option.m = pde.m;
option.rho = pde.rho;
option.gravity = pde.gravity;
option.eps_reg = getconfig(foInvConfig,'eps_reg',1e-3);
option.maxIt = 200;
option.tol = 1e-11;
option.residual_tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.edgequadorder = 4;
option.assemble_tangent = true;

%% Periodic P1 beta parameterization
Nm = Nx;
xBeta = (0:Nm-1)'*L/Nm;
xi = mod(xBeta,L)/L;

if 0
    betaTrue = 2*ones(size(xi));
    betaTrueName = 'constant';
elseif 0
    betaTrue = 2*(0.85+0.30*xi);
    betaTrueName = 'linear';
elseif 0
    betaTrue = 2*(0.9+0.4*(2*xi-1).^2);
    betaTrueName = 'quadratic';
elseif 1
    betaTrue = 2*(1+0.25*cos(2*pi*xi));
    betaTrueName = 'trigonometric';
elseif 0
    betaTrue = 2*(1+0.20*cos(2*pi*xi)+0.10*sin(4*pi*xi));
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
    betaPerturbation = 0.15*sin(2*pi*xi)+...
        0.08*cos(4*pi*xi)+0.05;
    perturbationName = 'mixed trigonometric';
end

betaInitial = betaTrue.*(1+betaPerturbation);
qTrue = log(betaTrue);
q = log(betaInitial);
betaPlotX = linspace(0,L,401)';

D = spdiags([-ones(Nm,1),ones(Nm,1)],[0,1],Nm,Nm);
D(Nm,1) = 1;

%% Synthetic top observation
[uTrue,eqnTrue,trueInfo] = solveforward(qTrue,[],pde,option,xBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');

dofNode = eqnTrue.dofNode;
surfaceLevel = -slope*dofNode(:,1);
tolGeometry = 100*eps(max(1,max(abs(node(:)))));
topDof = find(abs(dofNode(:,2)-surfaceLevel)<tolGeometry & ...
              dofNode(:,1)<L-tolGeometry);
[~,order] = sort(dofNode(topDof,1));
topDof = topDof(order);
xObs = dofNode(topDof,1);
topWeight = boundaryweights(xObs,L,slope);

dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

fprintf(['  beta true: %s, initial perturbation: %s, ',...
    'beta parameters: %d, top observations: %d\n'],...
    betaTrueName,perturbationName,Nm,numel(topDof));

%% Inverse options
maxInverseIt = getconfig(foInvConfig,'maxInverseIt',50);
alpha = getconfig(foInvConfig,'alpha',1e-11); % iteration 14
% alpha = getconfig(foInvConfig,'alpha',1e-12); % iteration 23
lambda = getconfig(foInvConfig,'lambda',1e-8);
pcgTolerance = 1e-8;
pcgMaxIt = 50;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;
useLineSearch = getconfig(foInvConfig,'useLineSearch',defaultlinesearch());

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
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,xBeta,L);
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
    G = assembleparameterderivative(eqn,q,xBeta,L,Nm);

    observationGradient = zeros(size(u));
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    observationGradientMaster = eqn.periodicProjection'*observationGradient;
    adjoint = eqn.tangent'\(-observationGradientMaster);
    gradient = G'*adjoint+alpha*(D'*(D*q));

    betaCurrent = exp(q);
    history.objective(k) = objective;
    history.dataResidual(k) = sqrt(...
        (topWeight'*(residual.^2))/dataNormSquared);
    history.parameterError(k) = norm(betaCurrent-betaTrue)/norm(betaTrue);
    history.parameterErrorLinf(k) = norm(betaCurrent-betaTrue,inf);
    history.parameterErrorRelativeLinf(k) = ...
        history.parameterErrorLinf(k)/norm(betaTrue,inf);
    history.gradientNorm(k) = norm(gradient);
    history.picardSteps(k) = forwardInfo.itStep;

    if k == 1
        derivativeCheck = verifyderivatives(q,u,eqn,G,gradient,...
            dataObs,topWeight,dataNormSquared,pde,option,xBeta,L,topDof,...
            D,alpha);
        optimizationForwardSolves = optimizationForwardSolves+...
            derivativeCheck.forwardSolves;
    end

    if norm(gradient) <= gradientTolerance
        printiterationrow(k,objective,history.parameterErrorLinf(k),...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,NaN,NaN,0,'grad');
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,D,alpha,lambda);
    [step,flag,relativeResidual,pcgIt] = pcg(...
        hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:FOAdjointPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        printiterationrow(k,objective,history.parameterErrorLinf(k),...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,pcgIt,relativeResidual,0,'step');
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
            [uTrial,~,trialInfo] = solveforward(...
                qTrial,u,pde,option,xBeta,L);
            optimizationForwardSolves = optimizationForwardSolves+1;
            if trialInfo.converged
                trialResidual = uTrial(topDof)-dataObs;
                trialRegularization = D*qTrial;
                trialObjective = 0.5*(topWeight'*...
                    (trialResidual.^2))/dataNormSquared+...
                    0.5*alpha*(trialRegularization'*...
                    trialRegularization);
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
            warning('iFEM:FOAdjointNoStep',...
                'No decreasing step was found; increasing LM damping.');
        end
    else
        qTrial = q+step;
        [uTrial,~,trialInfo] = solveforward(...
            qTrial,u,pde,option,xBeta,L);
        optimizationForwardSolves = optimizationForwardSolves+1;
        if trialInfo.converged
            q = qTrial;
            uWarm = uTrial;
            lambda = max(lambda/3,1e-12);
        else
            lambda = 10*lambda;
            warning('iFEM:FOAdjointNoStep',...
                'Undamped Gauss-Newton trial did not converge.');
        end
    end

    printiterationrow(k,objective,history.parameterErrorLinf(k),...
        history.parameterErrorRelativeLinf(k),norm(gradient),...
        forwardInfo.itStep,pcgIt,relativeResidual,lineSearchCount,'');

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = exp(q);
betaErrorLinf = norm(betaRecovered-betaTrue,inf);
betaErrorRelativeLinf = betaErrorLinf/norm(betaTrue,inf);
fprintf('\nSummary\n');
fprintf('  optimization forward solves: %d\n',optimizationForwardSolves);
fprintf('  final beta Linf error: %.04e, relative Linf: %.04e\n',...
    betaErrorLinf,betaErrorRelativeLinf);
fprintf(['  derivative check: state %.04e, grad %.04e, GN %.04e ',...
         '(FD %.04e, adj %.04e)\n'],...
    derivativeCheck.stateError,derivativeCheck.gradientError,...
    derivativeCheck.gaussNewtonError,derivativeCheck.finiteDifference,...
    derivativeCheck.adjointDirection);

figure(1);
betaTruePlot = periodicP1(betaPlotX,xBeta,betaTrue,L);
betaInitialPlot = periodicP1(betaPlotX,xBeta,betaInitial,L);
betaRecoveredPlot = periodicP1(betaPlotX,xBeta,betaRecovered,L);
plot(betaPlotX,betaTruePlot,'k-',...
    'LineWidth',1.8,'DisplayName','true');
hold on;
plot(betaPlotX,betaInitialPlot,'b--',...
    'LineWidth',1.2,'DisplayName','initial');
plot(betaPlotX,betaRecoveredPlot,'r-',...
    'LineWidth',1.4,'DisplayName','recovered');
plot(xBeta,betaRecovered,'ro','MarkerSize',5,...
    'DisplayName','recovered nodes');
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title(sprintf('FO sin-bed beta inversion: true %s, perturbation %s',...
    betaTrueName,perturbationName));

figure(2);
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
iteration = 1:numel(history.objective);
semilogy(iteration,history.objective,'o-',...
    'LineWidth',1.4,'DisplayName','objective');
grid on;
xlabel('inverse iteration');
ylabel('objective');
legend('Location','best');
title('objective history');

figure(4);
[uRecovered,~,~,~,~,~] = solveforward(q,uWarm,pde,option,xBeta,L);
plot(xObs,dataObs,'ko',xObs,uRecovered(topDof),'r-',...
    'LineWidth',1.4);
grid on;
xlabel('x');
ylabel('surface horizontal velocity');
legend('observation','recovered prediction','Location','best');

figure(5);
nNode = size(node,1);
subplot(1,2,1);
showmesh(node,elem);
title('FO sinusoidal-bed mesh','FontSize',14);

subplot(1,2,2);
trisurf(elem,node(:,1),node(:,2),uRecovered(1:nNode),...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('recovered u','FontSize',14);
view(2);
drawnow;

function [u,eqn,info,node,elem,bdFlag] = solveforward(...
        q,u0,pde,option,xBeta,L)
    beta = exp(q(:));
    option.beta = @(pt) periodicP1(pt(:,1),xBeta,beta,L);
    option.A = pde.A;
    option.n = pde.n;
    option.m = pde.m;
    option.rho = pde.rho;
    option.gravity = pde.gravity;
    option.f = pde.f;
    if isempty(u0)
        if isfield(option,'u0')
            option = rmfield(option,'u0');
        end
    else
        option.u0 = u0;
    end
    [soln,eqn,info,node,elem,bdFlag] = FirstOrderP2(option);
    u = soln.u;
end

function G = assembleparameterderivative(eqn,q,xBeta,L,Nm)
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        direction = zeros(Nm,1);
        direction(j) = 1;
        deltaBeta = beta.*direction;
        directionFunction = @(pt) periodicP1(pt(:,1),xBeta,deltaBeta,L);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,D,alpha,lambda)
    incrementalMaster = eqn.tangent\(-G*direction);
    incrementalState = eqn.periodicProjection*incrementalMaster;
    incrementalObservation = zeros(size(incrementalState));
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;
    incrementalObservationMaster = ...
        eqn.periodicProjection'*incrementalObservation;
    incrementalAdjoint = eqn.tangent'\(-incrementalObservationMaster);
    product = G'*incrementalAdjoint+alpha*(D'*(D*direction))+...
        lambda*direction;
end

function check = verifyderivatives(q,u,eqn,G,gradient,dataObs,topWeight,...
        dataNormSquared,pde,option,xBeta,L,topDof,D,alpha)
    direction = sin((1:numel(q))');
    direction = direction/norm(direction);
    epsilon = 1e-3;

    [uPlus,~,plusInfo] = solveforward(...
        q+epsilon*direction,u,pde,option,xBeta,L);
    [uMinus,~,minusInfo] = solveforward(...
        q-epsilon*direction,u,pde,option,xBeta,L);
    assert(plusInfo.converged && minusInfo.converged,...
        'Gradient-check forward solve failed.');

    plusResidual = uPlus(topDof)-dataObs;
    minusResidual = uMinus(topDof)-dataObs;
    plusRegularization = D*(q+epsilon*direction);
    minusRegularization = D*(q-epsilon*direction);
    plusObjective = 0.5*(topWeight'*(plusResidual.^2))/...
        dataNormSquared+0.5*alpha*...
        (plusRegularization'*plusRegularization);
    minusObjective = 0.5*(topWeight'*(minusResidual.^2))/...
        dataNormSquared+0.5*alpha*...
        (minusRegularization'*minusRegularization);
    finiteDifference = (plusObjective-minusObjective)/(2*epsilon);
    adjointDirection = gradient'*direction;
    relativeGradientError = abs(finiteDifference-adjointDirection)/...
        max([eps,abs(finiteDifference),abs(adjointDirection)]);

    incrementalMaster = eqn.tangent\(-G*direction);
    incrementalState = eqn.periodicProjection*incrementalMaster;
    finiteDifferenceState = (uPlus-uMinus)/(2*epsilon);
    relativeStateError = norm(finiteDifferenceState-incrementalState)/...
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

function printiterationheader()
    width = [3,12,12,11,12,7,5,10,2,4];
    label = {'it','objective','betaLinfAbs','betaLinfRel',...
        '|grad|','fPicard','pcgIt','pcgRel','ls','stop'};
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

function printiterationrow(k,objective,betaLinfAbs,betaLinfRel,...
        gradientNorm,forwardPicard,pcgIt,pcgRel,...
        lineSearchCount,stopReason)
    width = [3,12,12,11,12,7,5,10,2,4];
    if isnan(pcgIt), pcgItText = '-'; else, pcgItText = sprintf('%d',pcgIt); end
    if isnan(pcgRel), pcgRelText = '-'; else, pcgRelText = sprintf('%.04e',pcgRel); end
    if isempty(stopReason), stopReason = '-'; end
    value = {sprintf('%d',k),sprintf('%.04e',objective),...
        sprintf('%.04e',betaLinfAbs),sprintf('%.04e',betaLinfRel),...
        sprintf('%.04e',gradientNorm),sprintf('%d',forwardPicard),...
        pcgItText,pcgRelText,sprintf('%d',lineSearchCount),stopReason};
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

function history = trimhistory(history,k)
    fields = fieldnames(history);
    for j = 1:numel(fields)
        history.(fields{j}) = history.(fields{j})(1:k,:);
    end
end

function value = defaultlinesearch()
    value = false;
end

function value = getconfig(config,name,defaultValue)
    if isfield(config,name)
        value = config.(name);
    else
        value = defaultValue;
    end
end
