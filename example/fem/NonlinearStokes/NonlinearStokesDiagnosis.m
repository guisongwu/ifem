%% NONLINEARSTOKESDIAGNOSIS Fair n=1 versus n=3 comparison.
%
% Keep the geometry, mesh, bed parameter, initial guess, observations and
% optimization settings fixed.  Change only Glen's exponent n.

close all;
clear variables;

L = 1;
H = 0.5;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==0.5','Robin','y==0');
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

Nm = round(L/h);
assert(abs(Nm*h-L) <= 100*eps(max(1,L)),...
    'The parameter grid requires L/h to be an integer.');
xBeta = (0:Nm-1)'*h;
betaTrue = 1+0.1*cos(2*pi*xBeta/L);
betaInitial = betaTrue+0.1*(sin(2*pi*xBeta/L)+0.25);
qTrue = log(betaTrue);
qInitial = log(betaInitial);

pdeBase = struct;
pdeBase.A = 1;
pdeBase.m = 1/3;
pdeBase.rho = 1;
pdeBase.gravity = [0,-1];
pdeBase.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.assemble_tangent = false;

inverseOption.maxIt = 10;
inverseOption.fdStep = 1e-3;
inverseOption.lambda = 1e-7;
inverseOption.stepTolerance = 1e-7;

nList = [1,3];
result = repmat(struct,2,1);

for experiment = 1:numel(nList)
    pde = pdeBase;
    pde.n = nList(experiment);
    fprintf('\n========== Glen n = %g ==========\n',pde.n);

    [uTrue,trueInfo] = solveforward(qTrue,[],pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(trueInfo.converged,'Truth solve failed for n=%g.',pde.n);
    dataObs = uTrue(topDof);
    dataNormSquared = max(topWeight'*(dataObs.^2),eps);

    [uInitial,initialInfo] = solveforward(qInitial,[],pde,option,...
        node,elem,bdFlag,xBeta,L);
    assert(initialInfo.converged,'Initial solve failed for n=%g.',pde.n);
    [JInitial,fdInfo] = finiteDifferenceJacobian(qInitial,uInitial,...
        pde,option,node,elem,bdFlag,xBeta,L,topDof,topWeight,...
        dataNormSquared,inverseOption.fdStep);
    assert(all(fdInfo),'A perturbed solve failed for n=%g.',pde.n);

    [~,S,V] = svd(JInitial,'econ');
    singularValues = diag(S);
    parameterError = qInitial-qTrue;
    modeCoefficient = V'*parameterError;
    observableContribution = singularValues.*modeCoefficient;

    fprintf('initial relative beta error: %.6e\n',...
        norm(betaInitial-betaTrue)/norm(betaTrue));
    fprintf('initial normalized data residual: %.6e\n',...
        sqrt((topWeight'*((uInitial(topDof)-dataObs).^2))/...
        dataNormSquared));
    fprintf('singular values:\n');
    fprintf('  %.12e\n',singularValues);
    fprintf('condition(J): %.6e\n',...
        singularValues(1)/singularValues(end));
    fprintf('q-error coefficients in right singular vectors:\n');
    fprintf('  %.12e\n',modeCoefficient);
    fprintf('data contribution sigma_i*(v_i'' error):\n');
    fprintf('  %.12e\n',observableContribution);

    [qRecovered,history,JFinal] = invertparameter(qInitial,qTrue,dataObs,...
        topWeight,dataNormSquared,pde,option,inverseOption,node,elem,...
        bdFlag,xBeta,L,topDof);

    result(experiment).n = pde.n;
    result(experiment).singularValuesInitial = singularValues;
    result(experiment).conditionInitial = singularValues(1)/...
        singularValues(end);
    result(experiment).modeCoefficient = modeCoefficient;
    result(experiment).observableContribution = observableContribution;
    result(experiment).qRecovered = qRecovered;
    result(experiment).betaRecovered = exp(qRecovered);
    result(experiment).history = history;
    result(experiment).singularValuesFinal = svd(JFinal);

    fprintf('final relative beta error: %.6e\n',...
        history.parameterError(end));
    fprintf('final objective: %.6e\n',history.objective(end));
end

figure(1);
stairs([xBeta;L],[betaTrue;betaTrue(1)],'k-',...
    'LineWidth',1.8,'DisplayName','true');
hold on;
stairs([xBeta;L],[betaInitial;betaInitial(1)],'k--',...
    'LineWidth',1.2,'DisplayName','initial');
for experiment = 1:numel(result)
    beta = result(experiment).betaRecovered;
    stairs([xBeta;L],[beta;beta(1)],'-o','LineWidth',1.3,...
        'DisplayName',sprintf('recovered, n=%g',result(experiment).n));
end
hold off;
grid on;
xlabel('x');
ylabel('\beta');
legend('Location','best');
title('Fair comparison of bed-parameter recovery');

figure(2);
for experiment = 1:numel(result)
    semilogy(result(experiment).history.parameterError,'-o',...
        'LineWidth',1.3,...
        'DisplayName',sprintf('n=%g',result(experiment).n));
    hold on;
end
hold off;
grid on;
xlabel('inverse iteration');
ylabel('relative beta error');
legend('Location','best');

function [q,history,J] = invertparameter(q,qTrue,dataObs,topWeight,...
        dataNormSquared,pde,option,inverseOption,node,elem,bdFlag,...
        xBeta,L,topDof)
    Nm = numel(q);
    lambda = inverseOption.lambda;
    uWarm = [];

    history.objective = NaN(inverseOption.maxIt,1);
    history.parameterError = NaN(inverseOption.maxIt,1);
    history.dataResidual = NaN(inverseOption.maxIt,1);
    for k = 1:inverseOption.maxIt
        [u,info] = solveforward(q,uWarm,pde,option,...
            node,elem,bdFlag,xBeta,L);
        assert(info.converged,'Forward solve failed in iteration %d.',k);
        uWarm = u;
        residual = sqrt(topWeight/dataNormSquared).*...
            (u(topDof)-dataObs);
        objective = 0.5*(residual'*residual);

        history.objective(k) = objective;
        history.dataResidual(k) = norm(residual);
        history.parameterError(k) = norm(exp(q)-exp(qTrue))/...
            norm(exp(qTrue));

        [J,fdConverged] = finiteDifferenceJacobian(q,u,pde,option,...
            node,elem,bdFlag,xBeta,L,topDof,topWeight,...
            dataNormSquared,inverseOption.fdStep);
        assert(all(fdConverged),'A finite-difference solve failed.');

        gradient = J'*residual;
        normalMatrix = J'*J+lambda*speye(Nm);
        step = -normalMatrix\gradient;

        fprintf(['iteration %2d: objective %.6e, data residual %.6e, ',...
                 'beta error %.6e, lambda %.3e\n'],...
            k,objective,norm(residual),history.parameterError(k),lambda);

        if norm(step) <= inverseOption.stepTolerance*max(1,norm(q))
            history.objective = history.objective(1:k);
            history.dataResidual = history.dataResidual(1:k);
            history.parameterError = history.parameterError(1:k);
            return
        end

        accepted = false;
        stepLength = 1;
        for lineSearchIt = 1:10
            qTrial = q+stepLength*step;
            [uTrial,trialInfo] = solveforward(qTrial,u,pde,option,...
                node,elem,bdFlag,xBeta,L);
            if trialInfo.converged
                trialResidual = sqrt(topWeight/dataNormSquared).*...
                    (uTrial(topDof)-dataObs);
                trialObjective = 0.5*(trialResidual'*trialResidual);
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
        end
    end
end

function [J,converged] = finiteDifferenceJacobian(q,u,pde,option,...
        node,elem,bdFlag,xBeta,L,topDof,topWeight,dataNormSquared,...
        fdStep)
    Nm = numel(q);
    J = zeros(numel(topDof),Nm);
    converged = false(Nm,1);
    prediction = u(topDof);
    observationWeight = sqrt(topWeight/dataNormSquared);
    for j = 1:Nm
        qPerturbed = q;
        qPerturbed(j) = qPerturbed(j)+fdStep;
        [uPerturbed,info] = solveforward(qPerturbed,u,pde,option,...
            node,elem,bdFlag,xBeta,L);
        converged(j) = info.converged;
        J(:,j) = observationWeight.*...
            (uPerturbed(topDof)-prediction)/fdStep;
    end
end

function [u,info] = solveforward(q,u0,pde,option,...
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
