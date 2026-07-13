%% NSCONVERRATE Manufactured-solution convergence-rate test.

close all;
clear variables;

L = 1;
W = 1;
H = 0.5;
slope = 0.1;
eps_reg = 1e-2;
% hlist = [1/2;1/4;1/8];
hlist = [1/16];
nlevel = length(hlist);

errUx = zeros(nlevel,1);
errUy = zeros(nlevel,1);
errUz = zeros(nlevel,1);
errU = zeros(nlevel,1);
errP = zeros(nlevel,1);
iteration = zeros(nlevel,1);

pde = NSMMSData(L,W,H,slope,eps_reg);

option.periodic = true;
option.periodic_x = [0,L];
option.periodic_y = [0,W];
option.periodic_slope = [slope,0];
option.pressure_constraint = 'none';
option.eps_reg = eps_reg;
option.maxIt = 150;
option.tol = 1e-8;
option.residual_tol = 1e-8;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 5;
option.facequadorder = 5;

for level = 1:nlevel
    h = hlist(level);
    [node,elem] = cubemesh([0,1,0,1,0,1],h);
    bdFlag = setboundary3(node,elem,'Neumann','z==1',...
        'Dirichlet','z==0');

    node(:,1) = L*node(:,1);
    node(:,2) = W*node(:,2);
    zeta = node(:,3);
    node(:,3) = -slope*node(:,1)+H*zeta;

    [soln,~,info] = NonlinearStokes3P2P1(...
        node,elem,bdFlag,pde,option);

    errUx(level) = getL2error3(node,elem,pde.exactux,soln.ux,5);
    errUy(level) = getL2error3(node,elem,pde.exactuy,soln.uy,5);
    errUz(level) = getL2error3(node,elem,pde.exactuz,soln.uz,5);
    errU(level) = sqrt(errUx(level)^2+errUy(level)^2+...
        errUz(level)^2);
    pShift = pressureoffset(node,elem,soln.p,pde.exactp,5);
    errP(level) = getL2error3(node,elem,pde.exactp,soln.p-pShift,5);
    iteration(level) = info.itStep;

    fprintf('h=%7.5f  ||u-uh||=%9.3e  ||p-ph||=%9.3e  it=%d\n',...
        h,errU(level),errP(level),info.itStep);
end

rateU = [NaN;log2(errU(1:end-1)./errU(2:end))];
rateP = [NaN;log2(errP(1:end-1)./errP(2:end))];

fprintf('\n');
printsummary(hlist,errU,rateU,errP,rateP,iteration);

figure;
loglog(hlist,errU,'o-',hlist,errP,'s-','LineWidth',1.5);
set(gca,'XDir','reverse');
grid on;
xlabel('h');
ylabel('L^2 error');
legend('velocity','pressure','Location','northwest');
title('3-D manufactured-solution convergence');

function offset = pressureoffset(node,elem,ph,exactp,quadOrder)
[lambda,weight] = quadpts3(quadOrder);
[~,volume] = gradbasis3(node,elem);
integral = 0;
for q = 1:size(lambda,1)
    xq = lambda(q,1)*node(elem(:,1),:) + ...
         lambda(q,2)*node(elem(:,2),:) + ...
         lambda(q,3)*node(elem(:,3),:) + ...
         lambda(q,4)*node(elem(:,4),:);
    phq = lambda(q,1)*ph(elem(:,1)) + ...
          lambda(q,2)*ph(elem(:,2)) + ...
          lambda(q,3)*ph(elem(:,3)) + ...
          lambda(q,4)*ph(elem(:,4));
    integral = integral+sum(weight(q)*volume.*(phq-exactp(xq)));
end
offset = integral/sum(volume);
end

function printsummary(hlist,errU,rateU,errP,rateP,iteration)
fprintf('%10s  %12s  %7s  %12s  %7s  %4s\n',...
    'h','velocityL2','uRate','pressureL2','pRate','it');
for k = 1:length(hlist)
    fprintf('%10.5f  %12.4e  %7s  %12.4e  %7s  %4d\n',...
        hlist(k),errU(k),ratestr(rateU(k)),errP(k),...
        ratestr(rateP(k)),iteration(k));
end
end

function value = ratestr(rate)
if isnan(rate)
    value = '  --  ';
else
    value = sprintf('%7.2f',rate);
end
end
