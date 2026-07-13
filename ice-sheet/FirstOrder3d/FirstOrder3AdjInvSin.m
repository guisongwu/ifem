%% FIRSTORDER3ADJINVSIN 3-D first-order adjoint beta inversion on a sinusoidal bed.
%
% Synthetic top-surface horizontal velocity observations (u,v) are generated
% with the 3-D FO solver and used to recover q = log(beta) on the basal
% boundary.

close all;
clear variables;
set(groot,'DefaultFigureVisible','on');

%% Geometry and mesh
L = 5;
W = 5;
H = 1;
slope = 0.1;
bedAmplitude = 0.1*H;
Nx = 10;
Ny = 10;
Nz = 2;

fprintf(['3-D FO sinusoidal-bed beta inversion: L = %.04e, ',...
    'W = %.04e, H = %.04e\n'],L,W,H);

[refnode,elem] = cubemesh([0,1,0,1,0,1],[1/Nx,1/Ny,1/Nz]);
bdFlag = setboundary3(refnode,elem,'Neumann','z==1','Robin','z==0');
node = maptosinbed(refnode,L,W,H,slope,bedAmplitude);

[~,edge] = dof3P2(elem);
N = size(node,1);
Ndof = N+size(edge,1);
dofNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
tolGeometry = 1000*eps(max(1,max(abs(node(:)))));
topScalarDof = find(abs(dofNode(:,3)+slope*dofNode(:,1)-H) < ...
    tolGeometry & dofNode(:,1) < L-tolGeometry & ...
    dofNode(:,2) < W-tolGeometry);
topDof = [topScalarDof;Ndof+topScalarDof];
topWeightScalar = (L*W/numel(topScalarDof))*ones(numel(topScalarDof),1);
topWeight = [topWeightScalar;topWeightScalar];

%% FO model
pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = 1;
pde.gradS = @(pt) repmat([-slope,0],size(pt,1),1);
pde.g_N = [];

option.periodic = [1 2];
option.periodicBox = [0,L;0,W;min(node(:,3)),max(node(:,3))];
option.periodicSlope = [slope,0];
option.bed_condition = 'sliding';
option.eps_reg = 1e-3;
option.maxIt = 200;
option.tol = 1e-6;
option.residual_tol = 1e-6;
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
D = periodicgradient2d(NxBeta,NyBeta);

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
alpha = 0;
lambda = 1e-6;
pcgTolerance = 1e-7;
pcgMaxIt = 50;
gradientTolerance = 1e-9;
stepTolerance = 1e-8;
objectiveTolerance = 1e-15;
useLineSearch = defaultlinesearch();
uWarm = [];
optimizationForwardSolves = 0;

history.objective = NaN(maxInverseIt,1);
history.dataResidual = NaN(maxInverseIt,1);
history.parameterError = NaN(maxInverseIt,1);
history.parameterErrorLinf = NaN(maxInverseIt,1);
history.parameterErrorRelativeLinf = NaN(maxInverseIt,1);
history.gradientNorm = NaN(maxInverseIt,1);
history.picardSteps = NaN(maxInverseIt,1);
history.iterationTime = NaN(maxInverseIt,1);

fprintf(['\n it    objective      betaLinfRel   |grad|    ',...
    'fPicard pcgIt   pcgRel   ls    time(s)  stop\n']);
fprintf([' --  ------------  ------------  ------------  ------- ----- ',...
    '---------- --  --------  ----\n']);
for k = 1:maxInverseIt
    iterationStart = tic;
    [u,eqn,forwardInfo] = solveforward(q,uWarm,pde,option,...
        node,elem,bdFlag,xBeta,yBeta,L,W);
    optimizationForwardSolves = optimizationForwardSolves+1;
    assert(forwardInfo.converged,...
        'Forward solve failed at inverse iteration %d.',k);
    uWarm = u;

    residual = u(topDof)-dataObs;
    dataObjective = 0.5*(topWeight'*(residual.^2))/dataNormSquared;
    regularization = D*q;
    objective = dataObjective+0.5*alpha*(regularization'*regularization);
    G = assembleparameterderivative(eqn,q,xBeta,yBeta,L,W);

    observationGradient = zeros(size(u));
    observationGradient(topDof) = topWeight.*residual/dataNormSquared;
    observationGradientMaster = eqn.periodicProjection'*...
        observationGradient;
    adjoint = eqn.tangent'\(-observationGradientMaster);
    gradient = G'*adjoint+alpha*(D'*(D*q));

    betaCurrent = exp(q);
    betaTrueVector = exp(qTrue);
    history.objective(k) = objective;
    history.dataResidual(k) = sqrt((topWeight'*(residual.^2))/...
        dataNormSquared);
    history.parameterError(k) = norm(betaCurrent-betaTrueVector)/...
        norm(betaTrueVector);
    history.parameterErrorLinf(k) = norm(betaCurrent-betaTrueVector,inf);
    history.parameterErrorRelativeLinf(k) = ...
        history.parameterErrorLinf(k)/norm(betaTrueVector,inf);
    history.gradientNorm(k) = norm(gradient);
    history.picardSteps(k) = forwardInfo.itStep;

    if objective < objectiveTolerance
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,NaN,NaN,0,history.iterationTime(k),'obj');
        history = trimhistory(history,k);
        break
    end

    if norm(gradient) <= gradientTolerance
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,NaN,NaN,0,history.iterationTime(k),'grad');
        history = trimhistory(history,k);
        break
    end

    hessian = @(direction) gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,D,alpha,lambda);
    [step,flag,pcgRel,pcgIt] = pcg(...
        hessian,-gradient,pcgTolerance,pcgMaxIt);
    if flag ~= 0
        warning('iFEM:FO3AdjointPCG',...
            'PCG did not reach the requested tolerance.');
    end

    if norm(step) <= stepTolerance*max(1,norm(q))
        history.iterationTime(k) = toc(iterationStart);
        printiteration(k,objective,...
            history.parameterErrorRelativeLinf(k),norm(gradient),...
            forwardInfo.itStep,pcgIt,pcgRel,0,...
            history.iterationTime(k),'step');
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
            warning('iFEM:FO3AdjointNoStep',...
                'No decreasing step was found; increasing LM damping.');
        end
    else
        qTrial = q+step;
        [uTrial,~,trialInfo] = solveforward(qTrial,u,pde,option,...
            node,elem,bdFlag,xBeta,yBeta,L,W);
        optimizationForwardSolves = optimizationForwardSolves+1;
        if trialInfo.converged
            q = qTrial;
            uWarm = uTrial;
            lambda = max(lambda/3,1e-12);
        else
            lambda = 10*lambda;
            warning('iFEM:FO3AdjointNoStep',...
                'Undamped Gauss-Newton trial did not converge.');
        end
    end

    history.iterationTime(k) = toc(iterationStart);
    printiteration(k,objective,history.parameterErrorRelativeLinf(k),...
        norm(gradient),forwardInfo.itStep,pcgIt,pcgRel,lineSearchCount,...
        history.iterationTime(k),'-');

    if k == maxInverseIt
        history = trimhistory(history,k);
    end
end

betaRecovered = reshape(exp(q),NxBeta,NyBeta);
betaErrorLinf = norm(betaRecovered(:)-betaTrue(:),inf);
betaErrorRelativeLinf = betaErrorLinf/norm(betaTrue(:),inf);
fprintf('\nSummary\n');
fprintf('  optimization forward solves: %d\n',optimizationForwardSolves);
fprintf('  final beta Linf error: %.04e, relative Linf: %.04e\n',...
    betaErrorLinf,betaErrorRelativeLinf);

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
sgtitle('3-D FO sinusoidal-bed beta inversion');

figure(2);
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
semilogy(iteration,history.parameterErrorLinf,'s-','LineWidth',1.4,...
    'DisplayName','absolute beta Linf');
hold on;
semilogy(iteration,history.parameterErrorRelativeLinf,'^-',...
    'LineWidth',1.4,'DisplayName','relative beta Linf');
hold off;
grid on;
xlabel('inverse iteration');
ylabel('\beta error');
legend('Location','best');
title('\beta error history');

figure(3);
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

[uRecovered,~,~,solnRecovered] = solveforward(q,uWarm,pde,option,...
    node,elem,bdFlag,xBeta,yBeta,L,W);
figure(4);
set(gcf,'Visible','on');
clf;
plotsurfacefields(node,solnRecovered,Ndof,H,slope,L,W);
sgtitle('recovered top-surface FO fields');
exportepsfigures(mfilename);


function node = maptosinbed(refnode,L,W,H,slope,bedAmplitude)
    node = refnode;
    x = L*refnode(:,1);
    y = W*refnode(:,2);
    zeta = refnode(:,3);
    surface = H-slope*x;
    bed = -slope*x+bedAmplitude*sin(2*pi*x/L);
    node(:,1) = x;
    node(:,2) = y;
    node(:,3) = bed+zeta.*(surface-bed);
end

function plotsurfacefields(node,soln,Ndof,H,slope,L,W)
    top = abs(node(:,3)+slope*node(:,1)-H) < ...
        1000*eps(max(1,max(abs(node(:)))));
    topNode = node(top,:);
    topIndex = find(top);
    topU = soln.U(topIndex);
    topV = soln.U(Ndof+topIndex);
    topSpeed = sqrt(topU.^2+topV.^2);
    topElem = delaunay(topNode(:,1),topNode(:,2));
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile;
    plottopfield(topNode,topElem,topU,L,W,'surface u');
    nexttile;
    plottopfield(topNode,topElem,topV,L,W,'surface v');
    nexttile;
    plottopfield(topNode,topElem,topSpeed,L,W,'surface speed');
end

function plottopfield(topNode,topElem,value,L,W,titleText)
    trisurf(topElem,topNode(:,1),topNode(:,2),topNode(:,3),value,...
        'FaceColor','interp','EdgeColor','interp');
    axis equal tight;
    xlim([0,L]);
    ylim([0,W]);
    colorbar;
    xlabel('x');
    ylabel('y');
    title(titleText);
    view(2);
end

function product = gaussnewtonproduct(direction,eqn,G,topDof,...
        topWeight,dataNormSquared,D,alpha,lambda)
    incrementalMaster = eqn.tangent\(-G*direction);
    incrementalState = eqn.periodicProjection*incrementalMaster;
    incrementalObservation = zeros(size(incrementalState));
    incrementalObservation(topDof) = ...
        topWeight.*incrementalState(topDof)/dataNormSquared;
    incrementalObservationMaster = eqn.periodicProjection'*...
        incrementalObservation;
    incrementalAdjoint = eqn.tangent'\(-incrementalObservationMaster);
    product = G'*incrementalAdjoint+alpha*(D'*(D*direction))+...
        lambda*direction;
end

function G = assembleparameterderivative(eqn,q,xBeta,yBeta,L,W)
    Nm = numel(q);
    G = zeros(size(eqn.tangent,1),Nm);
    beta = exp(q(:));
    applyBetaDerivative = eqn.applyBetaDerivative;
    for j = 1:Nm
        deltaBeta = zeros(Nm,1);
        deltaBeta(j) = beta(j);
        betaDirection = reshape(deltaBeta,numel(xBeta),numel(yBeta));
        directionFunction = @(pt) periodicRectP1(pt(:,1),pt(:,2),...
            xBeta,yBeta,betaDirection,L,W);
        G(:,j) = applyBetaDerivative(directionFunction);
    end
end

function [u,eqn,info,soln] = solveforward(q,u0,pde,option,...
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
    [soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);
    u = soln.U;
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

function D = periodicgradient2d(nx,ny)
    nm = nx*ny;
    rows = [];
    cols = [];
    vals = [];
    cursor = 0;
    for iy = 1:ny
        for ix = 1:nx
            id = sub2ind([nx,ny],ix,iy);
            ixp = mod(ix,nx)+1;
            iyp = mod(iy,ny)+1;
            cursor = cursor+1;
            rows = [rows;cursor;cursor]; %#ok<AGROW>
            cols = [cols;id;sub2ind([nx,ny],ixp,iy)]; %#ok<AGROW>
            vals = [vals;-1;1]; %#ok<AGROW>
            cursor = cursor+1;
            rows = [rows;cursor;cursor]; %#ok<AGROW>
            cols = [cols;id;sub2ind([nx,ny],ix,iyp)]; %#ok<AGROW>
            vals = [vals;-1;1]; %#ok<AGROW>
        end
    end
    D = sparse(rows,cols,vals,2*nm,nm);
end

function printiteration(k,objective,betaRel,gradientNorm,...
        forwardPicard,pcgIt,pcgRel,lineSearchCount,iterationTime,stopReason)
    if isnan(pcgIt)
        pcgText = sprintf('%3s','-');
    else
        pcgText = sprintf('%3d',pcgIt);
    end
    if isnan(pcgRel)
        pcgRelText = sprintf('%10s','-');
    else
        pcgRelText = sprintf('%10.04e',pcgRel);
    end
    if isempty(stopReason)
        stopReason = '-';
    end
    fprintf(['%3d  %.04e    %.04e    %.04e    %5d   %s  %s  ',...
        '%2d  %8.2f  %s\n'],...
        k,objective,betaRel,gradientNorm,...
        forwardPicard,pcgText,pcgRelText,lineSearchCount,...
        iterationTime,stopReason);
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
