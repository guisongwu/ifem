%% NONLINEARSTOKESMMS Manufactured-solution convergence test.

close all;
clear variables;

L = 1;
H = 0.5;
slope = 0.1;
% A moderate regularization keeps the manufactured viscosity smooth.
% The physical demo can still use the smaller value 1e-4.
eps_reg = 1e-2;
hlist = [1/4;1/8;1/16;1/32];
nlevel = length(hlist);

errUx = zeros(nlevel,1);
errUz = zeros(nlevel,1);
errU = zeros(nlevel,1);
errP = zeros(nlevel,1);
iteration = zeros(nlevel,1);

pde = NonlinearStokesMMSData(L,H,slope,eps_reg);

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = eps_reg;
option.maxIt = 150;
option.tol = 1e-11;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 9;

for level = 1:nlevel
    h = hlist(level);
    [node,elem] = squaremesh([0,L,0,H],h);
    topBoundaryExpression = sprintf('y==%.17g',H);
    bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
        'Robin','y==0');
    node(:,2) = node(:,2)-slope*node(:,1);

    [soln,~,info] = NonlinearStokesP2P1(...
        node,elem,bdFlag,pde,option);

    errUx(level) = getL2error(node,elem,pde.exactux,soln.ux,6);
    errUz(level) = getL2error(node,elem,pde.exactuz,soln.uz,6);
    errU(level) = hypot(errUx(level),errUz(level));
    errP(level) = getL2error(node,elem,pde.exactp,soln.p,6);
    iteration(level) = info.itStep;

    fprintf('h=%7.5f  ||u-uh||=%9.3e  ||p-ph||=%9.3e  it=%d\n',...
        h,errU(level),errP(level),info.itStep);
end

rateU = [NaN;log(errU(1:end-1)./errU(2:end))/log(2)];
rateP = [NaN;log(errP(1:end-1)./errP(2:end))/log(2)];

fprintf('\n');
disp(table(hlist,errU,rateU,errP,rateP,iteration,...
    'VariableNames',{'h','velocityL2','velocityRate',...
                     'pressureL2','pressureRate','PicardSteps'}));

figure;
loglog(hlist,errU,'o-',hlist,errP,'s-','LineWidth',1.5);
set(gca,'XDir','reverse');
grid on;
xlabel('h');
ylabel('L^2 error');
legend('velocity','pressure','Location','northwest');
title('Manufactured-solution convergence');
