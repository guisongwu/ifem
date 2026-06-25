%% NONLINEARSTOKESBMMS Manufactured solution on an ISMIP-HOM B bed.
%
% This is a diagnostic construction for a sliding manufactured solution on
% a B-type bed.  With the current straight-edged triangular geometry it is
% not a clean convergence-rate test for a smooth sinusoidal bed.

close all;
clear variables;

warning('iFEM:NonlinearStokesBMMSDiagnosticOnly',...
    ['NonlinearStokesBMMS is a diagnostic B-bed sliding MMS. ',...
     'Do not use its table as a formal convergence-rate test.']);

L = 1;
H = 0.5;
bedAmplitude = 0.05;
slope = tan(0.5*pi/180);
eps_reg = 1e-2;
hlist = [1/4;1/8;1/16;1/32];
nlevel = length(hlist);

pde = NonlinearStokesBMMSData(L,H,bedAmplitude,slope,eps_reg);

errUx = zeros(nlevel,1);
errUz = zeros(nlevel,1);
errU = zeros(nlevel,1);
errP = zeros(nlevel,1);
iteration = zeros(nlevel,1);
converged = false(nlevel,1);

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = eps_reg;
option.maxIt = 200;
option.tol = 1e-10;
option.residual_tol = 1e-10;
option.damping = 0.8;
option.printlevel = 0;
option.quadorder = 6;

for level = 1:nlevel
    h = hlist(level);
    Nx = round(L/h);
    Nz = round(1/h);
    [node,elem] = rectanglemesh(L,1,Nx,Nz);
    bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
    node = maptoexperimentb(node,L,H,bedAmplitude,slope);

    [soln,~,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);

    errUx(level) = getL2error(node,elem,pde.exactux,soln.ux,6);
    errUz(level) = getL2error(node,elem,pde.exactuz,soln.uz,6);
    errU(level) = hypot(errUx(level),errUz(level));
    errP(level) = getL2error(node,elem,pde.exactp,soln.p,6);
    iteration(level) = info.itStep;
    converged(level) = info.converged;

    fprintf(['h=%7.5f  converged=%d  ||u-uh||=%9.3e  ',...
        '||p-ph||=%9.3e  it=%d\n'],...
        h,info.converged,errU(level),errP(level),info.itStep);
end

rateU = [NaN;log(errU(1:end-1)./errU(2:end))/log(2)];
rateP = [NaN;log(errP(1:end-1)./errP(2:end))/log(2)];

fprintf('\n');
disp(table(hlist,converged,errU,rateU,errP,rateP,iteration,...
    'VariableNames',{'h','converged','velocityL2','velocityRate',...
                     'pressureL2','pressureRate','PicardSteps'}));

figure(1);
loglog(hlist,errU,'o-',hlist,errP,'s-','LineWidth',1.5);
set(gca,'XDir','reverse');
grid on;
xlabel('h');
ylabel('L^2 error');
legend('velocity','pressure','Location','northwest');
title('B-bed sliding manufactured-solution convergence');

figure(2);
clf;
[node,elem] = rectanglemesh(L,1,round(L/hlist(end)),round(1/hlist(end)));
node = maptoexperimentb(node,L,H,bedAmplitude,slope);
[~,edge] = dofP2(elem);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
uxExact = pde.exactux(uNode);
uzExact = pde.exactuz(uNode);
subplot(2,2,1);
showmesh(node,elem);
title('B-type mesh');
subplot(2,2,2);
showresult(node,elem,uxExact(1:size(node,1)));
title('exact u_x');
subplot(2,2,3);
showresult(node,elem,uzExact(1:size(node,1)));
title('exact u_z');
subplot(2,2,4);
showresult(node,elem,pde.exactp(node));
title('exact p');

function node = maptoexperimentb(node,L,H,bedAmplitude,slope)
    x = node(:,1);
    r = node(:,2);
    surface = -slope*x;
    bed = surface-H+bedAmplitude*sin(2*pi*x/L);
    node(:,2) = bed+r.*(surface-bed);
end

function [node,elem] = rectanglemesh(L,H,Nx,Nz)
    x = linspace(0,L,Nx+1);
    z = linspace(0,H,Nz+1);
    [X,Z] = meshgrid(x,z);
    node = [X(:),Z(:)];

    cellId = reshape(1:(Nx+1)*(Nz+1),Nz+1,Nx+1);
    elem = zeros(2*Nx*Nz,3);
    cursor = 0;
    for ix = 1:Nx
        for iz = 1:Nz
            v1 = cellId(iz,ix);
            v2 = cellId(iz,ix+1);
            v3 = cellId(iz+1,ix);
            v4 = cellId(iz+1,ix+1);
            cursor = cursor+1;
            elem(cursor,:) = [v1,v2,v4];
            cursor = cursor+1;
            elem(cursor,:) = [v1,v4,v3];
        end
    end
end
