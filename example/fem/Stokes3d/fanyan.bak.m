%% NONLINEARSTOKES3ADJINVSLABBED 3-D adjoint beta inversion on a slab bed.
%
% This is a dimensionless 3-D counterpart of the 2-D slab-bed inverse
% example.  It recovers q = log(beta) from synthetic top-surface horizontal
% velocity observations (u,v) on a periodic sloping slab.

close all;
clear variables;
set(groot,'DefaultFigureVisible','on');

%% Geometry and mesh
L = 5;
W = 5;
H = 1;
slope = 0.1;
% Nx = 6;
% Ny = 6;
% Nz = 2;
Nx = 10;
Ny = 10;
Nz = 2;

fprintf('3-D slab-bed beta inversion: L = %.04e, W = %.04e, H = %.04e\n',...
    L,W,H);

[refnode,elem] = cubemesh([0,1,0,1,0,1],[1/Nx,1/Ny,1/Nz]);
bdFlag = setboundary3(refnode,elem,'Neumann','z==1','Robin','z==0');
node = maptoslab(refnode,L,W,H,slope);

[~,edge] = dof3P2(elem);
N = size(node,1);
Nu = N+size(edge,1);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
tolGeometry = 1000*eps(max(1,max(abs(node(:)))));
topScalarDof = find(abs(uNode(:,3)+slope*uNode(:,1)-H) < tolGeometry & ...
    uNode(:,1) < L-tolGeometry & uNode(:,2) < W-tolGeometry);
topDof = [topScalarDof;Nu+topScalarDof];
topWeightScalar = (L*W/numel(topScalarDof))*ones(numel(topScalarDof),1);
topWeight = [topWeightScalar;topWeightScalar];

%% Nonlinear Stokes model
pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = [0,0,-1];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.periodic_y = [0,W];
option.periodic_slope = [slope,0];
option.pressure_constraint = 'none';
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-10;
option.residual_tol = 1e-10;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 5;
option.facequadorder = 5;
option.assemble_tangent = true;

%% Periodic P1 beta parameterization
NxBeta = 10;
NyBeta = 10;
xBeta = (0:NxBeta-1)'*L/NxBeta;
yBeta = (0:NyBeta-1)'*W/NyBeta;
[XB,YB] = ndgrid(xBeta,yBeta);
xi = XB/L;
eta = YB/W;

betaTrue = 2*(1+0.20*cos(2*pi*xi).*cos(2*pi*eta));
betaInitial = betaTrue.*(1+0.15*sin(2*pi*xi).*cos(2*pi*eta)+0.05);
qTrue = log(betaTrue(:));
q = log(betaInitial(:));
Nm = numel(q);

fprintf('  beta parameters: %d, top velocity observations: %d\n',...
    Nm,numel(topDof));

%% Synthetic observation
[uTrue,~,trueInfo] = solveforward(qTrue,[],pde,option,...
    node,elem,bdFlag,xBeta,yBeta,L,W);
assert(trueInfo.converged,'The truth solve did not converge.');
dataObs = uTrue(topDof);
dataNormSquared = max(topWeight'*(dataObs.^2),eps);

%% Inversion loop
maxInverseIt = 20;
lambda = 1e-6;
pcgTolerance = 1e-8;
pcgMaxIt = 50;
gradientTolerance = 1e-9;
stepTolerance = 1e-8;
useLineSearch = defaultlinesearch();
uWarm = [];

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.parameterErrorLinf = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);
history.picardSteps = NaN(maxInverseIt,1);
history.iterationTime = NaN(maxInverseIt,1);

fprintf('\n it    objective      betaL2       betaLinf      |grad|    fPicard pcgIt ls    time(s)\n');
fprintf(' --  ------------  ------------  ------------  ------------  ------- ----- --  --------\n');
for k = 1:maxInverseIt
    iterationStart = tic;
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,yBeta,L,W);
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    objective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    G = assembleparameterderivative(eqn,q,xBeta,yBeta,L,W);

    observationGradient = zeros(size(eqn.tangent,1),1);
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    adjoint = eqn.tangent'\(-observationGradient);
    gradient = G'*adjoint;

    betaCurrent = exp(q);
    betaTrueVector = exp(qTrue);
    history.objective(k) = objective;
    history.dataResidual(k) = sqrt((topWeight'*(residual.^2))/...
        dataNormSquared);
    history.parameterError(k) = norm(betaCurrent-betaTrueVector)/...
        norm(betaTrueVector);
    history.parameterErrorLinf(k) = norm(betaCurrent-betaTrueVector,inf);
    history.gradientNorm(k) = norm(gradient);
    history.picardSteps(k) = forwardInfo.itStep;

    if objective < 1e-15
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,history.parameterError(k),...
            history.parameterErrorLinf(k),norm(gradient),...
            forwardInfo.itStep,NaN,0,history.iterationTime(k));
        history = trimhistory(history,k);
        break
    end

    if norm(gradient) <= gradientTolerance
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,history.parameterError(k),...
            history.parameterErrorLinf(k),norm(gradient),...
            forwardInfo.itStep,NaN,0,history.iterationTime(k));
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,...
        topDof,topWeight,dataNormSquared,lambda);
    [step,flag,~,pcgIt] = pcg(hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:NonlinearStokes3AdjInvPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,history.parameterError(k),...
            history.parameterErrorLinf(k),norm(gradient),...
            forwardInfo.itStep,pcgIt,0,history.iterationTime(k));
        history = trimhistory(history,k);
        break
    end

    lineSearchCount = 0;
    if useLineSearch
        accepted = false;
        stepLength = 1;
        for ls = 1:10
            lineSearchCount = lineSearchCount+1;
            qTrial = q+stepLength*step;
            [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
                node,elem,bdFlag,xBeta,yBeta,L,W);
            if trialInfo.converged
                trialResidual = uTrial(topDof)-dataObs;
                trialObjective = 0.5*(topWeight'*(trialResidual.^2))/...
                    dataNormSquared;
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
            warning('iFEM:NonlinearStokes3AdjInvNoStep',...
                'No decreasing step was found; increasing LM damping.');
        end
    else
        qTrial = q+step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,yBeta,L,W);
        if trialInfo.converged
            q = qTrial;
            uWarm = uTrial;
            lambda = max(lambda/3,1e-12);
        else
            lambda = 10*lambda;
            warning('iFEM:NonlinearStokes3AdjInvNoStep',...
                'Undamped Gauss-Newton trial did not converge.');
        end
    end

    history.iterationTime(k) = toc(iterationStart);
    printiteration(k,objective,history.parameterError(k),...
        history.parameterErrorLinf(k),norm(gradient),...
        forwardInfo.itStep,pcgIt,lineSearchCount,...
        history.iterationTime(k));

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = reshape(exp(q),NxBeta,NyBeta);
fprintf('\nFinal relative beta error: %.04e\n',...
    norm(betaRecovered(:)-betaTrue(:))/norm(betaTrue(:)));

figure(1);
set(gcf,'Visible','on');
clf;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
surf(XB,YB,betaTrue);
title('true beta');
xlabel('x'); ylabel('y'); zlabel('\beta');
nexttile;
surf(XB,YB,betaInitial);
title('initial beta');
xlabel('x'); ylabel('y'); zlabel('\beta');
nexttile;
surf(XB,YB,betaRecovered);
title('recovered beta');
xlabel('x'); ylabel('y'); zlabel('\beta');
sgtitle('3-D slab-bed beta inversion');

figure(2);
set(gcf,'Visible','on');
clf;
betaPlotX = linspace(0,L,401)';
ySliceIndex = [0,2,4,6,8]+1;
for kSlice = 1:numel(ySliceIndex)
    subplot(2,3,kSlice);
    iy = ySliceIndex(kSlice);
    betaTrueSlice = periodicRectP1(betaPlotX,...
        yBeta(iy)*ones(size(betaPlotX)),xBeta,yBeta,betaTrue,L,W);
    betaInitialSlice = periodicRectP1(betaPlotX,...
        yBeta(iy)*ones(size(betaPlotX)),xBeta,yBeta,betaInitial,L,W);
    betaRecoveredSlice = periodicRectP1(betaPlotX,...
        yBeta(iy)*ones(size(betaPlotX)),xBeta,yBeta,betaRecovered,L,W);
    plot(betaPlotX,betaTrueSlice,'k-','LineWidth',1.8,...
        'DisplayName','true');
    hold on;
    plot(betaPlotX,betaInitialSlice,'b--','LineWidth',1.2,...
        'DisplayName','initial');
    plot(betaPlotX,betaRecoveredSlice,'r-','LineWidth',1.4,...
        'DisplayName','recovered');
    plot(xBeta,betaRecovered(:,iy),'ro','MarkerSize',4,...
        'DisplayName','recovered nodes');
    hold off;
    grid on;
    xlabel('x');
    ylabel('\beta');
    title(sprintf('y = %.2f',yBeta(iy)));
    if kSlice == 1
        legend('Location','best');
    end
end
subplot(2,3,6);
axis off;
sgtitle('\beta slices');

figure(3);
set(gcf,'Visible','on');
clf;
iteration = 1:numel(history.objective);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
semilogy(iteration,history.objective,'o-','LineWidth',1.4,...
    'DisplayName','objective');
grid on;
xlabel('inverse iteration');
ylabel('objective');
legend('Location','best');
title('objective history');
nexttile;
semilogy(iteration,history.parameterError,'s-','LineWidth',1.4,...
    'DisplayName','relative beta L2');
hold on;
semilogy(iteration,history.parameterErrorLinf,'^-','LineWidth',1.4,...
    'DisplayName','absolute beta Linf');
hold off;
grid on;
xlabel('inverse iteration');
ylabel('\beta error');
legend('Location','best');
title('\beta error history');

[uRecovered,~,~,pRecovered] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,yBeta,L,W);
figure(4);
set(gcf,'Visible','on');
clf;
plotsurfacefields(node,uNode,uRecovered,pRecovered,Nu,H,slope,L,W);
sgtitle('recovered top-surface fields');

function node = maptoslab(refnode,L,W,H,slope)
    node = refnode;
    node(:,1) = L*refnode(:,1);
    node(:,2) = W*refnode(:,2);
    node(:,3) = H*refnode(:,3)-slope*node(:,1);
end

function plotsurfacefields(node,uNode,u,p,Nu,H,slope,L,W)
    topVelocity = abs(uNode(:,3)+slope*uNode(:,1)-H) < ...
        1000*eps(max(1,max(abs(uNode(:)))));
    topPressure = abs(node(:,3)+slope*node(:,1)-H) < ...
        1000*eps(max(1,max(abs(node(:)))));
    topVelocityNode = uNode(topVelocity,:);
    topPressureNode = node(topPressure,:);
    topU = u(topVelocity);
    topV = u(Nu+find(topVelocity));
    topP = p(topPressure);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile;
    plottopfield(topVelocityNode,topU,L,W,'surface u');
    nexttile;
    plottopfield(topVelocityNode,topV,L,W,'surface v');
    nexttile;
    plottopfield(topPressureNode,topP,L,W,'surface p');
end

function plottopfield(topNode,value,L,W,titleText)
    scatter(topNode(:,1),topNode(:,2),24,value,'filled');
    axis equal tight;
    xlim([0,L]);
    ylim([0,W]);
    colorbar;
    xlabel('x');
    ylabel('y');
    title(titleText);
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

function G = assembleparameterderivative(eqn,q,xBeta,yBeta,L,W)
    Nm = numel(q);
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    for j = 1:Nm
        deltaBeta = zeros(Nm,1);
        deltaBeta(j) = beta(j);
        directionFunction = @(pt) periodicRectP1(pt(:,1),pt(:,2),...
            xBeta,yBeta,reshape(deltaBeta,numel(xBeta),numel(yBeta)),...
            L,W);
        G(:,j) = eqn.applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,p] = solveforward(q,u0,pde,option,...
        node,elem,bdFlag,xBeta,yBeta,L,W)
    beta = reshape(exp(q(:)),numel(xBeta),numel(yBeta));
    pde.beta = @(pt) periodicRectP1(pt(:,1),pt(:,2),xBeta,yBeta,...
        beta,L,W);
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

function value = periodicRectP1(x,y,xNode,yNode,nodalValue,L,W)
    nx = numel(xNode);
    ny = numel(yNode);
    hx = L/nx;
    hy = W/ny;
    xWrapped = mod(x,L);
    yWrapped = mod(y,W);
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

function printiteration(k,objective,betaL2,betaLinf,gradientNorm,...
        forwardPicard,pcgIt,lineSearchCount,iterationTime)
    if isnan(pcgIt)
        pcgText = '  -';
    else
        pcgText = sprintf('%3d',pcgIt);
    end
    fprintf('%3d  %.04e    %.04e    %.04e    %.04e    %5d   %s  %2d  %8.2f\n',...
        k,objective,betaL2,betaLinf,gradientNorm,...
        forwardPicard,pcgText,lineSearchCount,iterationTime);
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
