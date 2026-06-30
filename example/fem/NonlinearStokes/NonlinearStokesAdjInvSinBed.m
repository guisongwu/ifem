%% NONLINEARSTOKESADJINVSINBED Invert beta on a sinusoidal bed.
%
% This script follows NonlinearStokesAdjInvSlabBed.m, but uses an
% ISMIP-HOM-B-like sinusoidal bed and a sliding basal boundary.
% Synthetic top-surface velocity data are generated on the same discrete
% geometry and then used to recover the periodic basal friction coefficient
% beta.  The inversion variable is q=log(beta).

close all;
clear variables;
set(groot,'DefaultFigureVisible','on');

%% Geometry and mesh
% Fixed-thickness cases.  Change L to compare how the horizontal
% length L affects beta identifiability from top-surface velocity data.
if 0
    L = 1;
elseif 0
    L = 2;
elseif 1
    L = 4;
end
H = 1;
slope = tan(0.5*pi/180);
bedAmplitude = 0.1*H;
h = 0.1;
Nx = max(4,round(L/h));
Nz = max(3,round(H/h));

fprintf(['Sin-bed beta inversion case: L = %.04e, H = %.04e, ',...
    'h = %.04e, Nx = %d, Nz = %d\n'],L,H,h,Nx,Nz);

[node,elem] = rectanglemesh(L,1,Nx,Nz);
% Boundary flags are set on the reference strip [0,L]x[0,1] before the
% B-type geometry map.  The physical top is y=-slope*x after mapping, but
% its reference location is always y=1, independent of H.
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node = maptoexperimentb(node,L,H,bedAmplitude,slope);

[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = -slope*uNode(:,1);
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
pde.m = 1;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-11;
option.residual_tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
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
elseif 0
    betaPerturbation = 0.20*sin(2*pi*xi)+0.05;
    perturbationName = 'trigonometric';
elseif 1
    betaPerturbation = 0.15*sin(2*pi*xi)+0.08*cos(4*pi*xi)+0.05;
    perturbationName = 'mixed trigonometric';
end

betaInitial = betaTrue.*(1+betaPerturbation);
qTrue = log(betaTrue);
q = log(betaInitial);

fprintf(['  beta true: %s, initial perturbation: %s, ',...
    'beta parameters: %d, top observations: %d\n'],...
    betaTrueName,perturbationName,Nm,numel(topDof));

%% Synthetic surface observation
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
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    optimizationForwardSolves = optimizationForwardSolves+1;
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    G = assembleparameterderivative(eqn,q,xBeta,L,Nm);

    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

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
        derivativeCheck = verifyderivatives(...
            q,u,eqn,G,gradient,dataObs,topWeight,dataNormSquared,...
            pde,option,node,elem,bdFlag,xBeta,L,topDof);
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
        topDof,topWeight,dataNormSquared,lambda);
    [step,flag,relativeResidual,pcgIt] = pcg(...
        hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:NonlinearAdjointPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        printiterationrow(k,objective,history.parameterErrorLinf(k),...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,pcgIt,relativeResidual,0,'step');
        history = trimhistory(history,k);
        break
    end

    accepted = false;
    stepLength = 1;
    lineSearchCount = 0;
    for lineSearchIt = 1:10
        lineSearchCount = lineSearchCount+1;
        qTrial = q+stepLength*step;
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
    'LineWidth',1.4,'DisplayName','recovered');
plot(xBeta,betaRecovered,'ro','MarkerSize',5,...
    'DisplayName','recovered nodes');
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title(sprintf('Beta inversion: true %s, perturbation %s',...
    betaTrueName,perturbationName));

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
title('recovered pressure','FontSize',14);
drawnow;

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
        history.(fields{j}) = history.(fields{j})(1:k);
    end
end
