%% NSSLABCONTINUATION Continuation for an ice slab.
%
% Reduce eps_reg one decade at a time.  Each converged velocity is used as
% the initial guess for the next nonlinear solve on the same mesh.

close all;
clear variables;

L = 5;
H = 1;
slope = 0.1;
h = 1/10;
% betaCase = 'constant';
% betaCase = 'linear';
% betaCase = 'quadratic';
betaCase = 'sin';
% betaCase = 'cos';
% betaCase = 'bump';
epsList = 10.^(-1:-1:-4);
nStage = length(epsList);

[node,elem] = squaremesh([0,L,0,H],h);
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

pde = struct;
pde.A = 1;
pde.n = 3;
pde.beta = makebedbeta(betaCase,L);
pde.m = 1/3;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.maxIt = 150;
option.tol = 1e-8;
option.damping = 0.8;
option.printlevel = 0;

converged = false(nStage,1);
iteration = NaN(nStage,1);
etaMin = NaN(nStage,1);
etaMax = NaN(nStage,1);
nCompleted = 0;

for stage = 1:nStage
    eps_reg = epsList(stage);
    option.eps_reg = eps_reg;

    if stage > 1
        option.u0 = soln.u;
    end

    [stageSoln,~,info] = NonlinearStokesP2P1(...
        node,elem,bdFlag,pde,option);

    converged(stage) = info.converged;
    iteration(stage) = info.itStep;
    etaMin(stage) = info.viscosityRange(end,1);
    etaMax(stage) = info.viscosityRange(end,2);

    fprintf('eps_reg=%8.1e  converged=%d  it=%3d  eta=[%9.3e,%9.3e]\n',...
        eps_reg,info.converged,info.itStep,etaMin(stage),etaMax(stage));

    if ~info.converged
        warning('iFEM:NSContinuationFailed',...
            'Continuation stopped at eps_reg = %.1e.',eps_reg);
        break
    end

    soln = stageSoln;
    nCompleted = stage;
end

result = table(epsList(:),converged,iteration,etaMin,etaMax,...
    'VariableNames',{'epsReg','converged','PicardSteps','etaMin','etaMax'});
disp(result(1:max(1,nCompleted+(nCompleted<nStage)),:));

if nCompleted > 0
    figure(1);
    set(gcf,'Visible','on');
    clf;
    subplot(2,2,1);
    showmesh(node,elem);
    axis equal;
    axis tight;
    title('mesh');

    subplot(2,2,2);
    trisurf(elem,node(:,1),node(:,2),soln.ux(1:size(node,1)),...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title(sprintf('horizontal velocity u_x, \\epsilon_{reg}=%.0e',...
        epsList(nCompleted)));
    view(2);

    subplot(2,2,3);
    trisurf(elem,node(:,1),node(:,2),soln.uz(1:size(node,1)),...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title(sprintf('vertical velocity u_z, \\epsilon_{reg}=%.0e',...
        epsList(nCompleted)));
    view(2);

    subplot(2,2,4);
    trisurf(elem,node(:,1),node(:,2),soln.p,...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title(sprintf('pressure, \\epsilon_{reg}=%.0e',...
        epsList(nCompleted)));
    view(2);

    figure(2);
    xPlot = linspace(0,L,401)';
    bedPoint = [xPlot,-slope*xPlot];
    betaPlot = evaluatebeta(pde.beta,bedPoint);
    plot(xPlot,betaPlot,'k-','LineWidth',1.6);
    grid on;
    xlabel('x');
    ylabel('\beta');
    title(sprintf('bed beta: %s',betaCase));
end

function beta = makebedbeta(betaCase,L)
    switch lower(betaCase)
        case 'constant'
            beta = 10;
        case 'linear'
            beta = @(pt) 8+4*periodicx(pt(:,1),L)/L;
        case 'quadratic'
            x = @(pt) periodicx(pt(:,1),L)/L;
            beta = @(pt) 8+8*(x(pt)-0.5).^2;
        case 'sin'
            beta = @(pt) 10+2*sin(2*pi*periodicx(pt(:,1),L)/L);
        case 'cos'
            beta = @(pt) 10+2*cos(2*pi*periodicx(pt(:,1),L)/L);
        case 'bump'
            x = @(pt) periodicx(pt(:,1),L)/L;
            beta = @(pt) 8+5*exp(-((x(pt)-0.5)/0.18).^2);
        otherwise
            error('Unknown betaCase: %s.',betaCase);
    end
end

function xWrapped = periodicx(x,L)
    xWrapped = mod(x,L);
end

function value = evaluatebeta(beta,pt)
    if isa(beta,'function_handle')
        value = beta(pt);
    else
        value = beta+zeros(size(pt,1),1);
    end
end
