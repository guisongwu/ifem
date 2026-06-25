%% NONLINEARSTOKESUNREGULARIZEDMMS Compare eps_reg against unregularized MMS.

close all;
clear variables;

L = 1;
H = .5;
slope = 0.1;
h = 1/16;
epsList = [1e-2;1e-4;1e-6;0];
nCase = length(epsList);

[node,elem] = squaremesh([0,L,0,H],h);
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);
pde = NonlinearStokesUnregularizedMMSData(L,H,slope);
[~,edge] = dofP2(elem);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
optionInitial = [pde.exactux(uNode);pde.exactuz(uNode)];

errUx = zeros(nCase,1);
errUz = zeros(nCase,1);
errP = zeros(nCase,1);
iteration = zeros(nCase,1);
converged = false(nCase,1);

option.periodic = true;
option.periodic_x = [0,L];
option.maxIt = 200;
option.tol = 1e-11;
option.residual_tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;
option.u0 = optionInitial;

for k = 1:nCase
    option.eps_reg = epsList(k);
    [soln,~,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);

    errUx(k) = getL2error(node,elem,pde.exactux,soln.ux,6);
    errUz(k) = getL2error(node,elem,pde.exactuz,soln.uz,6);
    errP(k) = getL2error(node,elem,pde.exactp,soln.p,6);
    iteration(k) = info.itStep;
    converged(k) = info.converged;

    fprintf(['eps_reg=%8.1e  converged=%d  it=%3d  ',...
        'errUx=%.4e  errUz=%.4e  errP=%.4e\n'],...
        epsList(k),info.converged,info.itStep,...
        errUx(k),errUz(k),errP(k));
end

result = table(epsList,converged,iteration,errUx,errUz,errP,...
    'VariableNames',{'epsReg','converged','PicardSteps',...
                     'uxL2','uzL2','pressureL2'});
disp(result);

figure;
loglog(max(epsList,realmin),errUx,'o-',...
    max(epsList,realmin),errUz,'s-',...
    max(epsList,realmin),errP,'^-','LineWidth',1.5);
grid on;
xlabel('\epsilon_{reg}');
ylabel('L^2 error');
legend('u_x','u_z','p','Location','best');
title('Errors against unregularized manufactured solution');
