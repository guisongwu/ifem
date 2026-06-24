%% NONLINEARSTOKESMMSCONTINUATION Regularization continuation for MMS.
%
% Reduce eps_reg one decade at a time.  Each converged velocity is used as
% the initial guess for the next nonlinear solve on the same mesh.

close all;
clear variables;

L = 1;
H = 1;
slope = 0.1;
h = 1/16;
epsList = 10.^(-1:-1:-4);
nStage = length(epsList);

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

option.periodic = true;
option.periodic_x = [0,L];
option.maxIt = 150;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 9;

converged = false(nStage,1);
iteration = NaN(nStage,1);
etaMin = NaN(nStage,1);
etaMax = NaN(nStage,1);
errU = NaN(nStage,1);
errP = NaN(nStage,1);
nCompleted = 0;

for stage = 1:nStage
    eps_reg = epsList(stage);
    pde = NonlinearStokesMMSData(L,H,slope,eps_reg);
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
    errUx = getL2error(node,elem,pde.exactux,stageSoln.ux,6);
    errUz = getL2error(node,elem,pde.exactuz,stageSoln.uz,6);
    errU(stage) = hypot(errUx,errUz);
    errP(stage) = getL2error(node,elem,pde.exactp,stageSoln.p,6);

    fprintf(['eps_reg=%8.1e  converged=%d  it=%3d  ',...
             'eta=[%9.3e,%9.3e]  ||u-uh||=%9.3e  ||p-ph||=%9.3e\n'],...
        eps_reg,info.converged,info.itStep,etaMin(stage),etaMax(stage),...
        errU(stage),errP(stage));

    if ~info.converged
        warning('iFEM:NonlinearStokesContinuationFailed',...
            'Continuation stopped at eps_reg = %.1e.',eps_reg);
        break
    end

    soln = stageSoln;
    nCompleted = stage;
end

result = table(epsList(:),converged,iteration,etaMin,etaMax,errU,errP,...
    'VariableNames',{'epsReg','converged','PicardSteps',...
                     'etaMin','etaMax','velocityL2','pressureL2'});
disp(result(1:max(1,nCompleted+(nCompleted<nStage)),:));

if nCompleted > 0
    figure(1);
    showresult(node,elem,soln.p);
    title(sprintf('MMS pressure, \\epsilon_{reg} = %.0e',...
        epsList(nCompleted)));

    figure(2);
    showresult(node,elem,soln.ux(1:size(node,1)));
    title(sprintf('MMS horizontal velocity, \\epsilon_{reg} = %.0e',...
        epsList(nCompleted)));
end
