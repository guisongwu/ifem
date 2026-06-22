%% NONLINEAR_STOKES_BETA_INVERSION Invert a spatially varying bed parameter.
%
% Recover the positive Weertman coefficient beta(x) from synthetic surface
% velocity observations.  Each forward evaluation solves the full nonlinear
% Stokes problem with Glen viscosity and a nonlinear sliding law.
%
% The inversion variable is q = log(beta), which guarantees beta > 0.
% A finite-difference Jacobian and a regularized Gauss-Newton/Levenberg-
% Marquardt step are used.  This is deliberately independent of the frozen
% Picard matrix returned by NonlinearStokesP2P1_periodic: that matrix is not
% the consistent Jacobian required by a nonlinear adjoint calculation.

close all;
clear variables;

%% Geometry and finite-element mesh
L = 1;
H = 1;
slope = 0.2;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

%% Nonlinear Stokes and Weertman-law data
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
option.maxIt = 150;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;

%% Periodic P1 representation of beta on the bed
% Keep the inverse space coarser than the state mesh.  Surface-only data
% do not reliably identify bed oscillations at the finite-element scale.
Nm = 4;
xBeta = (0:Nm-1)'*L/Nm;

betaTrue = 1 + 0.3*cos(2*pi*xBeta/L) ...
             - 0.1*sin(2*pi*xBeta/L);
qTrue = log(betaTrue);

betaInitial = ones(Nm,1);
q = log(betaInitial);

%% Surface observation degrees of freedom
[~,edge] = dofP2(elem);
N = size(node,1);
Nu = N + size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
surfaceLevel = H-slope*uNode(:,1);
tol = 100*eps(max(1,max(abs(node(:)))));
topDof = find(abs(uNode(:,2)-surfaceLevel)<tol ...
            & uNode(:,1)<L-tol);
[~,order] = sort(uNode(topDof,1));
topDof = topDof(order);
xObs = uNode(topDof,1);

% Use both surface velocity components.  Set observeVertical=false to use
% horizontal velocity only, as in the original linear inversion scripts.
observeVertical = true;

%% Generate synthetic observations
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,L);
assert(trueInfo.converged,...
    'The nonlinear solve for the synthetic truth did not converge.');
dataTrue = observation(uTrue,topDof,Nu,observeVertical);

noiseLevel = 0;
rng(1);
noiseScale = noiseLevel*max(norm(dataTrue)/sqrt(numel(dataTrue)),eps);
dataObs = dataTrue + noiseScale*randn(size(dataTrue));
dataScale = max(norm(dataObs)/sqrt(numel(dataObs)),eps);

%% Regularized Gauss-Newton/Levenberg-Marquardt inversion
maxInverseIt = 12;
fdStep = 1e-2;
lambda = 1e-6;
alpha = 1e-8;
relativeStepTolerance = 1e-5;

% Periodic first-difference regularization for q.
D = spdiags([-ones(Nm,1),ones(Nm,1)],[0,1],Nm,Nm);
D(Nm,1) = 1;
qReference = log(betaInitial);

history.objective = NaN(maxInverseIt,1);
history.dataMisfit = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);

uWarm = [];
for k = 1:maxInverseIt
    [u,~,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(forwardInfo.converged,...
        'The nonlinear forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    prediction = observation(u,topDof,Nu,observeVertical);
    residual = (prediction-dataObs)/dataScale;
    regularization = D*(q-qReference);
    dataMisfit = 0.5*(residual'*residual);
    objective = dataMisfit + 0.5*alpha*(regularization'*regularization);

    history.objective(k) = objective;
    history.dataMisfit(k) = dataMisfit;
    history.parameterError(k) = norm(exp(q)-betaTrue)/norm(betaTrue);

    fprintf(['inverse %2d: objective %.3e, data %.3e, ',...
             'relative beta error %.3e, forward steps %d\n'],...
        k,objective,dataMisfit,history.parameterError(k),...
        forwardInfo.itStep);

    % Forward finite differences of the complete nonlinear parameter-to-
    % observation map.  Warm starts reduce the Picard iteration count.
    J = zeros(numel(residual),Nm);
    for j = 1:Nm
        qPerturbed = q;
        qPerturbed(j) = qPerturbed(j)+fdStep;
        [uPerturbed,~,perturbedInfo] = solveforward(qPerturbed,u,pde,...
            option,node,elem,bdFlag,xBeta,L);
        assert(perturbedInfo.converged,...
            'Perturbed nonlinear solve failed for parameter %d.',j);
        perturbedData = observation(...
            uPerturbed,topDof,Nu,observeVertical);
        J(:,j) = (perturbedData-prediction)/(fdStep*dataScale);
    end

    gradient = J'*residual + alpha*(D'*(D*(q-qReference)));
    normalMatrix = J'*J + alpha*(D'*D) + lambda*speye(Nm);
    step = -normalMatrix\gradient;

    if norm(step) <= relativeStepTolerance*max(1,norm(q))
        history.objective = history.objective(1:k);
        history.dataMisfit = history.dataMisfit(1:k);
        history.parameterError = history.parameterError(1:k);
        break
    end

    % Backtracking accepts only a decrease of the full objective.
    accepted = false;
    stepLength = 1;
    for lineSearchIt = 1:8
        qTrial = q+stepLength*step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,L);
        if trialInfo.converged
            trialResidual = (observation(...
                uTrial,topDof,Nu,observeVertical)-dataObs)/dataScale;
            trialRegularization = D*(qTrial-qReference);
            trialObjective = 0.5*(trialResidual'*trialResidual) ...
                + 0.5*alpha*(trialRegularization'*trialRegularization);
            if trialObjective < objective
                q = qTrial;
                uWarm = uTrial;
                lambda = max(lambda/3,1e-10);
                accepted = true;
                break
            end
        end
        stepLength = stepLength/2;
    end

    if ~accepted
        lambda = 10*lambda;
        warning('iFEM:NonlinearStokesInversionNoStep',...
            'No decreasing step was found; increasing LM damping.');
        if lambda > 1e8
            history.objective = history.objective(1:k);
            history.dataMisfit = history.dataMisfit(1:k);
            history.parameterError = history.parameterError(1:k);
            break
        end
    end

    if k == maxInverseIt
        history.objective = history.objective(1:k);
        history.dataMisfit = history.dataMisfit(1:k);
        history.parameterError = history.parameterError(1:k);
    end
end

betaRecovered = exp(q);

%% Results
figure(1);
stairs([xBeta;L],[betaTrue;betaTrue(1)],'k-','LineWidth',1.8);
hold on;
stairs([xBeta;L],[betaInitial;betaInitial(1)],'b--','LineWidth',1.2);
stairs([xBeta;L],[betaRecovered;betaRecovered(1)],'r-o',...
    'LineWidth',1.4,'MarkerSize',5);
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('true','initial','recovered','Location','best');
title('Nonlinear Stokes bed-parameter inversion');

figure(2);
semilogy(1:numel(history.objective),history.objective,'o-',...
    1:numel(history.dataMisfit),history.dataMisfit,'s-','LineWidth',1.4);
grid on;
xlabel('inverse iteration');
ylabel('functional value');
legend('objective','data misfit','Location','best');

figure(3);
recoveredData = observation(uWarm,topDof,Nu,observeVertical);
plot(xObs,dataObs(1:numel(topDof)),'ko',...
     xObs,recoveredData(1:numel(topDof)),...
     'r-','LineWidth',1.4);
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
    [soln,~,info] = NonlinearStokesP2P1_periodic(...
        node,elem,bdFlag,pde,option);
    u = soln.u;
    p = soln.p;
end

function value = periodicP1(x,xNode,nodalValue,L)
    xWrapped = mod(x,L);
    value = interp1([xNode;L],[nodalValue;nodalValue(1)],...
        xWrapped,'linear');
end

function data = observation(u,topDof,Nu,observeVertical)
    ux = u(topDof);
    if observeVertical
        uz = u(Nu+topDof);
        data = [ux;uz];
    else
        data = ux;
    end
end
