%% NSFDINVERSION Finite-difference beta inversion on a slab bed.
%
% This script uses the same slab-bed test cases, boundary-integral
% objective, beta truth choices, initial perturbation choices, and
% visualization style as FullStokesAdjInvSlab.m.
% Its intentional algorithmic difference is that centered finite
% differences of the complete reduced objective provide both the gradient
% and full Hessian.

close all;
clear variables;
set(groot,'DefaultFigureVisible','on');

%% Geometry and mesh
% Slab cases.  H is fixed; toggle exactly one L branch to compare how
% horizontal length affects beta identifiability.
if 0
    L = 1;
elseif 0
    L = 2;
elseif 1
    L = 4;
end
H = 1;
h = 1/4;
slope = 0.1;

fprintf('FD slab-bed beta inversion case: L = %.04e, H = %.04e, h = %.04e\n',...
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
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 4;

%% Periodic P1 beta parameterization
Nm = round(L/h);
assert(abs(Nm*h-L) <= 100*eps(max(1,L)),...
    'The parameter grid requires L/h to be an integer.');
xBeta = (0:Nm-1)'*h;
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
maxInverseIt = 10;
alpha = 0;
fdStep = 1e-2;
lambda = 1e-7;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;
D = spdiags([-ones(Nm,1),ones(Nm,1)],[0,1],Nm,Nm);
D(Nm,1) = 1;

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.parameterErrorLinf = NaN(maxInverseIt,1);
history.parameterErrorRelativeLinf = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);
history.picardSteps = NaN(maxInverseIt,1);

uWarm = [];
optimizationForwardSolves = 0;
printiterationheader();
for k = 1:maxInverseIt
    [u,~,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    optimizationForwardSolves = optimizationForwardSolves+1;
    assert(forwardInfo.converged,...
        'The nonlinear forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    regularization = D*q;
    objective = dataObjective+0.5*alpha*(regularization'*regularization);

    betaCurrent = exp(q);
    history.objective(k) = objective;
    history.dataResidual(k) = sqrt(...
        (topWeight'*(residual.^2))/dataNormSquared);
    history.parameterError(k) = norm(betaCurrent-betaTrue)/norm(betaTrue);
    history.parameterErrorLinf(k) = norm(betaCurrent-betaTrue,inf);
    history.parameterErrorRelativeLinf(k) = ...
        history.parameterErrorLinf(k)/norm(betaTrue,inf);
    history.picardSteps(k) = forwardInfo.itStep;

    % Centered finite differences of the complete reduced objective.
    % The plus/minus values are reused by the gradient and Hessian diagonal.
    differenceStep = fdStep*max(1,abs(q));
    objectivePlus = zeros(Nm,1);
    objectiveMinus = zeros(Nm,1);
    gradient = zeros(Nm,1);
    hessian = zeros(Nm,Nm);
    derivativeForwardEvaluations = 0;
    derivativePicardSteps = 0;

    for i = 1:Nm
        qPlus = q;
        qMinus = q;
        qPlus(i) = qPlus(i)+differenceStep(i);
        qMinus(i) = qMinus(i)-differenceStep(i);

        [objectivePlus(i),~,plusInfo] = reducedobjective(...
            qPlus,u,pde,option,node,elem,bdFlag,xBeta,L,...
            topDof,dataObs,topWeight,dataNormSquared,D,alpha);
        [objectiveMinus(i),~,minusInfo] = reducedobjective(...
            qMinus,u,pde,option,node,elem,bdFlag,xBeta,L,...
            topDof,dataObs,topWeight,dataNormSquared,D,alpha);
        assert(plusInfo.converged && minusInfo.converged,...
            'Finite-difference solve failed for parameter %d.',i);

        derivativeForwardEvaluations = derivativeForwardEvaluations+2;
        derivativePicardSteps = derivativePicardSteps ...
            + plusInfo.itStep+minusInfo.itStep;
        gradient(i) = (objectivePlus(i)-objectiveMinus(i)) ...
            /(2*differenceStep(i));
        hessian(i,i) = (objectivePlus(i)-2*objective ...
            + objectiveMinus(i))/differenceStep(i)^2;
    end

    % Four-point centered stencil for each mixed Hessian entry.
    for i = 1:Nm-1
        for j = i+1:Nm
            mixedObjective = zeros(2,2);
            signs = [-1,1];
            for si = 1:2
                for sj = 1:2
                    qMixed = q;
                    qMixed(i) = qMixed(i) ...
                        + signs(si)*differenceStep(i);
                    qMixed(j) = qMixed(j) ...
                        + signs(sj)*differenceStep(j);
                    [mixedObjective(si,sj),~,mixedInfo] = ...
                        reducedobjective(...
                            qMixed,u,pde,option,node,elem,bdFlag,...
                            xBeta,L,topDof,dataObs,topWeight,...
                            dataNormSquared,D,alpha);
                    assert(mixedInfo.converged,...
                        ['Mixed finite-difference solve failed for ',...
                         'parameters %d and %d.'],i,j);
                    derivativeForwardEvaluations = ...
                        derivativeForwardEvaluations+1;
                    derivativePicardSteps = derivativePicardSteps ...
                        + mixedInfo.itStep;
                end
            end
            hessian(i,j) = (mixedObjective(2,2) ...
                - mixedObjective(2,1)-mixedObjective(1,2) ...
                + mixedObjective(1,1)) ...
                /(4*differenceStep(i)*differenceStep(j));
            hessian(j,i) = hessian(i,j);
        end
    end
    optimizationForwardSolves = optimizationForwardSolves ...
        + derivativeForwardEvaluations;

    rawSymmetryError = norm(hessian-hessian','fro') ...
        /max(1,norm(hessian,'fro'));
    hessian = 0.5*(hessian+hessian');
    hessianEigenvalues = sort(eig(hessian));
    largestEigenvalue = max(abs(hessianEigenvalues));
    smallestEigenvalue = min(abs(hessianEigenvalues));
    spectralCondition = largestEigenvalue ...
        /max(smallestEigenvalue,eps*max(1,largestEigenvalue));
    negativeEigenvalues = sum(hessianEigenvalues < ...
        -sqrt(eps)*max(1,largestEigenvalue));

    history.gradientNorm(k) = norm(gradient);
    if norm(gradient) <= gradientTolerance
        printiterationrow(k,objective,history.parameterErrorLinf(k),...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,derivativeForwardEvaluations,...
            derivativePicardSteps,0,'grad');
        history = trimhistory(history,k);
        break
    end

    step = -(hessian+lambda*speye(Nm))\gradient;

    if norm(step) <= stepTolerance*max(1,norm(q))
        printiterationrow(k,objective,history.parameterErrorLinf(k),...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,derivativeForwardEvaluations,...
            derivativePicardSteps,0,'step');
        history = trimhistory(history,k);
        break
    end

    accepted = false;
    stepLength = 1;
    lineSearchCount = 0;
    for lineSearchIt = 1:10
        lineSearchCount = lineSearchCount+1;
        qTrial = q+stepLength*step;
        [trialObjective,uTrial,trialInfo] = reducedobjective(...
            qTrial,u,pde,option,node,elem,bdFlag,xBeta,L,...
            topDof,dataObs,topWeight,dataNormSquared,D,alpha);
        optimizationForwardSolves = optimizationForwardSolves+1;
        if trialInfo.converged && trialObjective < objective
            q = qTrial;
            uWarm = uTrial;
            lambda = max(lambda/3,1e-12);
            accepted = true;
            break
        end
        stepLength = stepLength/2;
    end

    if ~accepted
        lambda = 10*lambda;
        warning('iFEM:NSInversionNoStep',...
            'No decreasing step was found; increasing LM damping.');
    end

    printiterationrow(k,objective,history.parameterErrorLinf(k),...
        history.parameterErrorRelativeLinf(k),norm(gradient),...
        forwardInfo.itStep,derivativeForwardEvaluations,...
        derivativePicardSteps,lineSearchCount,'');

    fprintf(['  FD Hessian symmetry %.04e, negative eigs %d, ',...
        'condition %.04e\n'],...
        rawSymmetryError,negativeEigenvalues,spectralCondition);

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

%% Results
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
    'LineWidth',1.4,'DisplayName','FD recovered');
plot(xBeta,betaRecovered,'ro','MarkerSize',5,...
    'DisplayName','recovered nodes');
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title(sprintf('FD slab-bed beta inversion: true %s, perturbation %s',...
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
[uRecovered,pRecovered] = solveforward(q,uWarm,pde,option,...
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

%% Local functions
exportepsfigures(mfilename);

function printiterationheader()
    width = [3,12,12,11,12,7,8,8,2,4];
    label = {'it','objective','betaLinfAbs','betaLinfRel',...
        '|grad|','fPicard','fdSolves','fdPicard','ls','stop'};
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
        gradientNorm,forwardPicard,fdSolves,fdPicard,...
        lineSearchCount,stopReason)
    width = [3,12,12,11,12,7,8,8,2,4];
    if isempty(stopReason)
        stopReason = '-';
    end
    value = {sprintf('%d',k),sprintf('%.04e',objective),...
        sprintf('%.04e',betaLinfAbs),sprintf('%.04e',betaLinfRel),...
        sprintf('%.04e',gradientNorm),sprintf('%d',forwardPicard),...
        sprintf('%d',fdSolves),sprintf('%d',fdPicard),...
        sprintf('%d',lineSearchCount),stopReason};
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

function [u,p,info] = solveforward(q,u0,pde,option,...
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
    [soln,~,info] = NonlinearStokesP2P1(...
        node,elem,bdFlag,pde,option);
    u = soln.u;
    p = soln.p;
end

function [objective,u,info] = reducedobjective(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L,topDof,dataObs,topWeight,...
        dataNormSquared,D,alpha)
    [u,~,info] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L);
    if ~info.converged
        objective = Inf;
        return
    end
    residual = u(topDof)-dataObs;
    regularization = D*q;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared+...
        0.5*alpha*(regularization'*regularization);
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
