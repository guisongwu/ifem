%% RECTANGULARNONLINEARSTOKES3ADJINV Adjoint beta inversion in a cuboid.
%
% Synthetic top-surface horizontal velocity data are generated with a known
% basal friction beta(x,y).  The inversion variable is q=log(beta) on a
% tensor-product parameter grid on the flat bed z=0.

close all;
clear variables;

%% Geometry and mesh
L = 2;
W = 1;
H = 1;
h = 0.5;
[node,elem] = cubemesh([0,L,0,W,0,H],h);
bdFlag = setboundary3(node,elem,'Dirichlet','all',...
    'Neumann',sprintf('z==%.17g',H),'Robin','z==0');

[~,edge] = dof3P2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
tolGeometry = 100*eps(max(1,max(abs(node(:)))));
topScalarDof = find(abs(uNode(:,3)-H)<tolGeometry);
topDof = topScalarDof;
topWeight = (L*W/numel(topDof))*ones(numel(topDof),1);

%% Nonlinear Stokes model
pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = [0.15,0,-1];
pde.g_N = @(pt) zeros(size(pt,1),3);

option.eps_reg = 1e-3;
option.maxIt = 80;
option.tol = 1e-9;
option.residual_tol = 1e-9;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 4;
option.facequadorder = 4;
option.assemble_tangent = true;

%% Bed beta parameterization
Nx = 5;
Ny = 4;
xBeta = linspace(0,L,Nx)';
yBeta = linspace(0,W,Ny)';
[XB,YB] = ndgrid(xBeta,yBeta);
betaTrue = 2*(1+0.25*cos(2*pi*XB/L).*cos(2*pi*YB/W));
betaInitial = betaTrue.*(1+0.15*sin(pi*XB/L).*cos(pi*YB/W));
qTrue = log(betaTrue(:));
q = log(betaInitial(:));
Nm = numel(q);

fprintf('3-D cuboid beta inversion: beta parameters %d, top observations %d\n',...
    Nm,numel(topDof));

%% Synthetic observation
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,yBeta);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

%% Inversion loop
maxInverseIt = 8;
lambda = 1e-6;
pcgTolerance = 1e-8;
pcgMaxIt = 40;
uWarm = [];

history.objective = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);

fprintf('\n it    objective     relBetaErr      |grad|    fPicard pcgIt\n');
fprintf(' --  ------------  ------------  ------------  ------- -----\n');
for k = 1:maxInverseIt
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,yBeta);
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    G = assembleparameterderivative(eqn,q,xBeta,yBeta);

    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

    betaCurrent = exp(q);
    history.objective(k) = objective;
    history.parameterError(k) = norm(betaCurrent-exp(qTrue))/norm(exp(qTrue));
    history.gradientNorm(k) = norm(gradient);

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    [step,flag,~,pcgIt] = pcg(hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:NonlinearStokes3AdjInvPCG',...
            'PCG did not reach the requested tolerance.');
    end

    accepted = false;
    stepLength = 1;
    for ls = 1:8
        qTrial = q+stepLength*step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,yBeta);
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
    end

    fprintf('%3d  %.04e    %.04e    %.04e    %5d   %3d\n',...
        k,objective,history.parameterError(k),norm(gradient),...
        forwardInfo.itStep,pcgIt);
    if norm(gradient) < 1e-9
        history = trimhistory(history,k);
        break
    end
    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = reshape(exp(q),Nx,Ny);
fprintf('\nFinal relative beta error: %.04e\n',...
    norm(betaRecovered(:)-betaTrue(:))/norm(betaTrue(:)));

figure(1);
subplot(1,3,1);
surf(XB,YB,betaTrue);
title('true beta');
xlabel('x'); ylabel('y'); zlabel('\beta');
subplot(1,3,2);
surf(XB,YB,betaInitial);
title('initial beta');
xlabel('x'); ylabel('y'); zlabel('\beta');
subplot(1,3,3);
surf(XB,YB,betaRecovered);
title('recovered beta');
xlabel('x'); ylabel('y'); zlabel('\beta');

figure(2);
semilogy(1:numel(history.objective),history.objective,'o-',...
    1:numel(history.parameterError),history.parameterError,'s-');
grid on;
legend('objective','relative beta error','Location','best');
xlabel('inverse iteration');

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,lambda)
    incrementalState = eqn.tangent\(-G*direction);
    incrementalObservation = zeros(size(eqn.tangent,1),1);
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;
    incrementalAdjoint = eqn.tangent'\(-incrementalObservation);
    product = G'*incrementalAdjoint+lambda*direction;
end

function G = assembleparameterderivative(eqn,q,xBeta,yBeta)
    Nm = numel(q);
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        deltaBeta = zeros(Nm,1);
        deltaBeta(j) = beta(j);
        directionFunction = @(pt) rectP1(pt(:,1),pt(:,2),...
            xBeta,yBeta,reshape(deltaBeta,numel(xBeta),numel(yBeta)));
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,yBeta)
    beta = reshape(exp(q(:)),numel(xBeta),numel(yBeta));
    pde.beta = @(pt) rectP1(pt(:,1),pt(:,2),xBeta,yBeta,beta);
    if isempty(u0)
        if isfield(option,'u0')
            option = rmfield(option,'u0');
        end
    else
        option.u0 = u0;
    end
    [soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);
    u = soln.u;
    p = soln.p;
end

function value = rectP1(x,y,xNode,yNode,nodalValue)
    value = interp2(yNode,xNode,nodalValue,y,x,'linear');
end

function history = trimhistory(history,k)
    fields = fieldnames(history);
    for j = 1:numel(fields)
        history.(fields{j}) = history.(fields{j})(1:k,:);
    end
end
