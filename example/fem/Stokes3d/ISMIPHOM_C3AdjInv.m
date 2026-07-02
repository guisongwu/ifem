%% ISMIPHOM_C3ADJINV Adjoint beta inversion with x/y periodic boundaries.
%
% This is an inverse-crime test on the 3-D ISMIP-HOM-C geometry.  Synthetic
% top-surface horizontal velocity data are generated from a positive
% periodic basal friction field beta(x,y), then q=log(beta) is recovered
% from top u_x observations.

close all;
clear variables;

%% Geometry and mesh
rho = 910;
gravity = 9.81;
alpha = 0.1*pi/180;
slope = tan(alpha);
L = 5000;
H = 1000;

Nx = 4;
Ny = 4;
Nz = 2;
[node,elem] = cubemesh([0,1,0,1,0,1],[1/Nx,1/Ny,1/Nz]);
bdFlag = setboundary3(node,elem,'Neumann','z==1','Robin','z==0');

node(:,1) = L*node(:,1);
node(:,2) = L*node(:,2);
zeta = node(:,3);
node(:,3) = -node(:,1)*slope-H+zeta*H;

[~,edge] = dof3P2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
tolGeometry = 1000*eps(max(1,max(abs(node(:)))));
topScalarDof = find(abs(uNode(:,3)+slope*uNode(:,1)) < tolGeometry & ...
    uNode(:,1) < L-tolGeometry & uNode(:,2) < L-tolGeometry);
topDof = topScalarDof;
topWeight = (L^2/numel(topDof))*ones(numel(topDof),1);

%% Nonlinear Stokes model
pde = struct;
pde.A = 1e-16;
pde.n = 3;
pde.m = 1;
pde.rho = rho;
pde.gravity = [0,0,-gravity];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.periodic_y = [0,L];
option.periodic_slope_x = slope;
option.pressure_constraint = 'none';
option.eps_reg = 1e-10;
option.maxIt = 120;
option.tol = 1e-8;
option.residual_tol = 1e-8;
option.damping = 0.7;
option.printlevel = 0;
option.quadorder = 4;
option.facequadorder = 4;
option.assemble_tangent = true;

%% Periodic beta parameterization
NxBeta = 4;
NyBeta = 4;
xBeta = (0:NxBeta-1)'*L/NxBeta;
yBeta = (0:NyBeta-1)'*L/NyBeta;
[XB,YB] = ndgrid(xBeta,yBeta);

% Positive C-like friction field.  Exact ISMIP-HOM C reaches beta=0, which
% cannot be represented by q=log(beta).
betaTrue = 1000+800*sin(2*pi*XB/L).*sin(2*pi*YB/L);
betaInitial = betaTrue.*(1+0.20*cos(2*pi*XB/L).*cos(2*pi*YB/L));
qTrue = log(betaTrue(:));
q = log(betaInitial(:));
Nm = numel(q);

fprintf(['3-D periodic C beta inversion: parameters %d, ',...
    'top observations %d\n'],Nm,numel(topDof));

%% Synthetic observation
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,yBeta,L);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

%% Inversion loop
maxInverseIt = 6;
lambda = 1e-4;
pcgTolerance = 1e-7;
pcgMaxIt = 30;
gradientTolerance = 1e-8;
uWarm = [];

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);

fprintf('\n it    objective      dataRel      betaRel       |grad|    fPicard pcgIt ls\n');
fprintf(' --  ------------  ------------  ------------  ------------  ------- ----- --\n');
for k = 1:maxInverseIt
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,yBeta,L);
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    G = assembleparameterderivative(eqn,q,xBeta,yBeta,L);

    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

    betaCurrent = exp(q);
    history.objective(k) = objective;
    history.dataResidual(k) = sqrt((topWeight'*(residual.^2))/dataNormSquared);
    history.parameterError(k) = norm(betaCurrent-exp(qTrue))/norm(exp(qTrue));
    history.gradientNorm(k) = norm(gradient);

    if norm(gradient) < gradientTolerance
        fprintf('%3d  %.04e    %.04e    %.04e    %.04e    %5d     -   0\n',...
            k,objective,history.dataResidual(k),history.parameterError(k),...
            norm(gradient),forwardInfo.itStep);
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    [step,flag,~,pcgIt] = pcg(hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:ISMIPHOMC3AdjInvPCG',...
            'PCG did not reach the requested tolerance.');
    end

    accepted = false;
    stepLength = 1;
    lineSearchCount = 0;
    for ls = 1:8
        lineSearchCount = lineSearchCount+1;
        qTrial = q+stepLength*step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,yBeta,L);
        if trialInfo.converged
            trialResidual = uTrial(topDof)-dataObs;
            trialObjective = 0.5*(topWeight'*(trialResidual.^2)) ...
                /dataNormSquared;
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
    end

    fprintf('%3d  %.04e    %.04e    %.04e    %.04e    %5d   %3d  %2d\n',...
        k,objective,history.dataResidual(k),history.parameterError(k),...
        norm(gradient),forwardInfo.itStep,pcgIt,lineSearchCount);

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = reshape(exp(q),NxBeta,NyBeta);
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

function G = assembleparameterderivative(eqn,q,xBeta,yBeta,L)
    Nm = numel(q);
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        deltaBeta = zeros(Nm,1);
        deltaBeta(j) = beta(j);
        directionFunction = @(pt) periodicRectP1(pt(:,1),pt(:,2),...
            xBeta,yBeta,reshape(deltaBeta,numel(xBeta),numel(yBeta)),L);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,yBeta,L)
    beta = reshape(exp(q(:)),numel(xBeta),numel(yBeta));
    pde.beta = @(pt) periodicRectP1(pt(:,1),pt(:,2),xBeta,yBeta,beta,L);
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

function value = periodicRectP1(x,y,xNode,yNode,nodalValue,L)
    nx = numel(xNode);
    ny = numel(yNode);
    hx = L/nx;
    hy = L/ny;
    xWrapped = mod(x,L);
    yWrapped = mod(y,L);
    ix0 = floor(xWrapped/hx)+1;
    iy0 = floor(yWrapped/hy)+1;
    tx = (xWrapped-xNode(ix0))/hx;
    ty = (yWrapped-yNode(iy0))/hy;
    ix1 = mod(ix0,nx)+1;
    iy1 = mod(iy0,ny)+1;
    v00 = nodalValue(sub2ind([nx,ny],ix0,iy0));
    v10 = nodalValue(sub2ind([nx,ny],ix1,iy0));
    v01 = nodalValue(sub2ind([nx,ny],ix0,iy1));
    v11 = nodalValue(sub2ind([nx,ny],ix1,iy1));
    value = (1-tx).*(1-ty).*v00+tx.*(1-ty).*v10+...
        (1-tx).*ty.*v01+tx.*ty.*v11;
end

function history = trimhistory(history,k)
    fields = fieldnames(history);
    for j = 1:numel(fields)
        history.(fields{j}) = history.(fields{j})(1:k,:);
    end
end
