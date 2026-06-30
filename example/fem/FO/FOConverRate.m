%% FOCONVERRATE Manufactured-solution convergence-rate test.

close all;
clear variables;

L = 4;
H = 1;
slope = 0.1;

% Use n = 1 for a clean P2 convergence benchmark.  The manufactured
% solution is written in the mapped vertical coordinate q = z+slope*x.
hlist = [1/2;1/4;1/8;1/16];
nlevel = length(hlist);

errU = zeros(nlevel,1);
errW = zeros(nlevel,1);
errP = zeros(nlevel,1);
iteration = zeros(nlevel,1);

pde = fommsdata(L,H,slope);

option.L = L;
option.H = H;
option.slope = slope;
option.periodic_x = [0,L];
option.A = pde.A;
option.n = pde.n;
option.m = pde.m;
option.beta = pde.beta;
option.rho = pde.rho;
option.gravity = pde.gravity;
option.eps_reg = pde.eps_reg;
option.f = pde.f;
option.maxIt = 20;
option.tol = 1e-12;
option.residual_tol = 1e-12;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 9;

fineSoln = [];
fineEqn = [];
fineNode = [];
fineElem = [];

for level = 1:nlevel
    h = hlist(level);
    option.h = [h,h];

    [soln,eqn,info,node,elem] = FirstOrderP2(option);

    errU(level) = getp2l2error(node,elem,eqn.elem2dof,pde.exactu,soln.u);
    errW(level) = getp2l2error(node,elem,eqn.elem2dof,pde.exactw,soln.w);
    errP(level) = getp2l2error(node,elem,eqn.elem2dof,pde.exactp,soln.p);
    iteration(level) = info.itStep;

    fprintf('h=%7.5f  ||u-uh||=%9.3e  ||w-wh||=%9.3e  ||p-ph||=%9.3e  it=%d\n',...
        h,errU(level),errW(level),errP(level),info.itStep);

    if level == nlevel
        fineSoln = soln;
        fineEqn = eqn;
        fineNode = node;
        fineElem = elem;
    end
end

rateU = [NaN;log(errU(1:end-1)./errU(2:end))/log(2)];
rateW = [NaN;log(errW(1:end-1)./errW(2:end))/log(2)];
rateP = [NaN;log(errP(1:end-1)./errP(2:end))/log(2)];

fprintf('\n');
printsummary(hlist,errU,rateU,errW,rateW,errP,rateP,iteration);

figure;
loglog(hlist,errU,'o-',hlist,errW,'d-',hlist,errP,'s-','LineWidth',1.5);
set(gca,'XDir','reverse');
grid on;
xlabel('h');
ylabel('L^2 error');
legend('u','w','pressure','Location','northwest');
title('FO manufactured-solution convergence');

plotdiagnostics(fineNode,fineElem,fineEqn,fineSoln,L,H,slope);

function pde = fommsdata(L,H,slope)
pde.A = 1;
pde.n = 1;
pde.m = 1;
pde.beta = 10;
pde.rho = 1;
pde.gravity = 1;
pde.eps_reg = 1e-10;
pde.f = @(pt) force(pt,L,H,slope);
pde.exactu = @(pt) exactu(pt,L,H,slope);
pde.exactw = @(pt) exactw(pt,L,H,slope);
pde.exactp = @(pt) exactp(pt,L,H,slope,pde.rho,pde.gravity);
end

function printsummary(hlist,errU,rateU,errW,rateW,errP,rateP,iteration)
fprintf('%10s  %12s  %7s  %12s  %7s  %12s  %7s  %4s\n',...
    'h','uL2','uRate','wL2','wRate','pL2','pRate','it');
for k = 1:length(hlist)
    fprintf('%10.5f  %12.4e  %7s  %12.4e  %7s  %12.4e  %7s  %4d\n',...
        hlist(k),errU(k),ratestr(rateU(k)),errW(k),ratestr(rateW(k)),...
        errP(k),ratestr(rateP(k)),iteration(k));
end
end

function value = ratestr(rate)
if isnan(rate)
    value = '  --  ';
else
    value = sprintf('%7.2f',rate);
end
end

function value = exactu(pt,L,H,slope)
x = pt(:,1);
q = pt(:,2)+slope*x;
k = 2*pi/L;
value = sin(k*x).*shape(q,H);
end

function value = exactw(pt,L,H,slope)
x = pt(:,1);
q = pt(:,2)+slope*x;
k = 2*pi/L;
[G,Gint] = shapeintegral(q,H);
value = -(k*cos(k*x).*Gint + slope*sin(k*x).*G);
end

function value = exactp(pt,L,H,slope,rho,gravity)
x = pt(:,1);
q = pt(:,2)+slope*x;
k = 2*pi/L;
[G,Gp] = shape(q,H);
ux = k*cos(k*x).*G + slope*sin(k*x).*Gp;
value = rho*gravity*(H-q) - ux;
end

function value = force(pt,L,H,slope)
x = pt(:,1);
q = pt(:,2)+slope*x;
k = 2*pi/L;
[G,Gp,Gpp] = shape(q,H);
uxx = -k^2*sin(k*x).*G + 2*slope*k*cos(k*x).*Gp + ...
      slope^2*sin(k*x).*Gpp;
uzz = sin(k*x).*Gpp;
value = -(2*uxx + 0.5*uzz);
end

function [G,Gp,Gpp] = shape(q,H)
G = q.^2.*(H-q).^2;
Gp = 2*q.*(H-q).^2 - 2*q.^2.*(H-q);
Gpp = 2*H^2 - 12*H*q + 12*q.^2;
end

function [G,Gint] = shapeintegral(q,H)
G = q.^2.*(H-q).^2;
Gint = H^2*q.^3/3 - H*q.^4/2 + q.^5/5;
end

function errL2 = getp2l2error(node,elem,elem2dof,exact,value)
[lambda,weight] = quadpts(9);
errL2sq = 0;
for q = 1:size(lambda,1)
    phi = p2basis2(lambda(q,:));
    xq = lambda(q,1)*node(elem(:,1),:) + ...
         lambda(q,2)*node(elem(:,2),:) + ...
         lambda(q,3)*node(elem(:,3),:);
    valueq = value(elem2dof)*phi';
    exactq = exact(xq);
    [~,area] = gradbasis(node,elem);
    errL2sq = errL2sq + sum(weight(q)*area.*(valueq-exactq).^2);
end
errL2 = sqrt(errL2sq);
end

function phi = p2basis2(l)
phi = [l(1)*(2*l(1)-1), l(2)*(2*l(2)-1), ...
       l(3)*(2*l(3)-1), 4*l(2)*l(3), ...
       4*l(3)*l(1), 4*l(1)*l(2)];
end

function plotdiagnostics(node,elem,eqn,soln,L,H,slope)
nNode = size(node,1);

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
showmesh(node,elem);
axis equal;
axis tight;
title('mesh');

nexttile;
trisurf(elem,node(:,1),node(:,2),soln.u(1:nNode),...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('horizontal velocity u');
view(2);

nexttile;
trisurf(elem,node(:,1),node(:,2),soln.w(1:nNode),...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('vertical velocity w');
view(2);

nexttile;
trisurf(elem,node(:,1),node(:,2),soln.p(1:nNode),...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('pressure p');
view(2);

sgtitle(sprintf('FO finest-grid fields, L = %.3g, H = %.3g, slope = %.3g',...
    L,H,slope));
end
