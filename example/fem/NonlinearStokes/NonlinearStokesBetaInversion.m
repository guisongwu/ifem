%% NONLINEARSTOKESBETAINVERSION Finite-difference inversion for bed beta(x).
%
% This script uses the same model, observations, initial parameter,
% objective, regularization, damping, line search, and stopping tolerances
% as NonlinearStokesAdjointInversion.  Its only intentional algorithmic
% difference is that centered finite differences of the complete reduced
% objective provide both the gradient and full Hessian.

close all;
clear variables;

%% Geometry and mesh
L = 1;
H = 0.5;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==0.5','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

%% Surface observation degrees of freedom
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

%% Periodic P1 parameterization
Nm = 4;
xBeta = (0:Nm-1)'*L/Nm;
betaTrue = 1+0.1*cos(2*pi*xBeta/L);
betaInitial = betaTrue+0.1*(sin(2*pi*xBeta/L)+0.25);
qTrue = log(betaTrue);
q = log(betaInitial);
qReference = q;

D = spdiags([-ones(Nm,1),ones(Nm,1)],[0,1],Nm,Nm);
D(Nm,1) = 1;

%% Synthetic surface observation
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataScale = max(norm(dataObs)/sqrt(numel(dataObs)),eps);

%% Inverse options
maxInverseIt = 10;
fdStep = 1e-2;
alpha = 1e-8;
lambda = 1e-7;
stepTolerance = 1e-7;
gradientTolerance = 1e-9;

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);

uWarm = [];
for k = 1:maxInverseIt
    [u,~,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(forwardInfo.converged,...
        'The nonlinear forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = (u(topDof)-dataObs)/dataScale;
    regularization = D*(q-qReference);
    objective = 0.5*(residual'*residual) ...
        + 0.5*alpha*(regularization'*regularization);

    history.objective(k) = objective;
    history.dataResidual(k) = norm(residual);
    history.parameterError(k) = norm(exp(q)-betaTrue)/norm(betaTrue);

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
            topDof,dataObs,dataScale,D,qReference,alpha);
        [objectiveMinus(i),~,minusInfo] = reducedobjective(...
            qMinus,u,pde,option,node,elem,bdFlag,xBeta,L,...
            topDof,dataObs,dataScale,D,qReference,alpha);
        assert(plusInfo.converged && minusInfo.converged,...
            'Finite-difference solve failed for parameter %d.',i);

        derivativeForwardEvaluations = ...
            derivativeForwardEvaluations+2;
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
                            xBeta,L,topDof,dataObs,dataScale,D,...
                            qReference,alpha);
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

    fprintf(['  FD gradient %.3e, Hessian symmetry %.3e, ',...
             'derivative forward solves %d, Picard steps %d\n'],...
        norm(gradient),rawSymmetryError,...
        derivativeForwardEvaluations,derivativePicardSteps);
    fprintf('  Hessian eigenvalues:');
    fprintf(' %.3e',hessianEigenvalues);
    fprintf('\n');
    fprintf(['  Hessian negative eigenvalues %d, ',...
             'spectral condition %.3e\n'],...
        negativeEigenvalues,spectralCondition);
    history.gradientNorm(k) = norm(gradient);
    fprintf(['FD inverse %2d: objective %.6e, data %.6e, ',...
             'beta error %.6e, |gradient| %.6e\n'],...
        k,objective,norm(residual),history.parameterError(k),...
        norm(gradient));

    if norm(gradient) <= gradientTolerance
        history = trimhistory(history,k);
        break
    end

    step = -(hessian+lambda*speye(Nm))\gradient;

    if norm(step) <= stepTolerance*max(1,norm(q))
        history = trimhistory(history,k);
        break
    end

    accepted = false;
    stepLength = 1;
    for lineSearchIt = 1:10
        qTrial = q+stepLength*step;
        [trialObjective,uTrial,trialInfo] = reducedobjective(...
            qTrial,u,pde,option,node,elem,bdFlag,xBeta,L,...
            topDof,dataObs,dataScale,D,qReference,alpha);
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
        warning('iFEM:NonlinearStokesInversionNoStep',...
            'No decreasing step was found; increasing LM damping.');
    end

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = exp(q);

%% Results
figure(1);
stairs([xBeta;L],[betaTrue;betaTrue(1)],'k-',...
    'LineWidth',1.8,'DisplayName','true');
hold on;
stairs([xBeta;L],[betaInitial;betaInitial(1)],'b--',...
    'LineWidth',1.2,'DisplayName','initial');
stairs([xBeta;L],[betaRecovered;betaRecovered(1)],'r-o',...
    'LineWidth',1.4,'DisplayName','finite-difference recovered');
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title('Nonlinear Stokes finite-difference inversion');

figure(2);
semilogy(history.parameterError,'o-','LineWidth',1.4);
grid on;
xlabel('inverse iteration');
ylabel('relative beta error');

figure(3);
[uRecovered,~,~] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,L);
plot(xObs,dataObs,'ko',xObs,uRecovered(topDof),'r-',...
    'LineWidth',1.4);
grid on;
xlabel('x');
ylabel('surface horizontal velocity');
legend('observation','recovered prediction','Location','best');

%% Local functions
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
        node,elem,bdFlag,xBeta,L,topDof,dataObs,dataScale,D,...
        qReference,alpha)
    [u,~,info] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,L);
    if ~info.converged
        objective = Inf;
        return
    end
    residual = (u(topDof)-dataObs)/dataScale;
    regularization = D*(q-qReference);
    objective = 0.5*(residual'*residual) ...
        + 0.5*alpha*(regularization'*regularization);
end

function value = periodicP1(x,xNode,nodalValue,L)
    xWrapped = mod(x,L);
    value = interp1([xNode;L],[nodalValue;nodalValue(1)],...
        xWrapped,'linear');
end

function history = trimhistory(history,k)
    fields = fieldnames(history);
    for j = 1:numel(fields)
        history.(fields{j}) = history.(fields{j})(1:k);
    end
end
