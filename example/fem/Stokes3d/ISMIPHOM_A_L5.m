%% ISMIPHOM_A_L5 3-D ISMIP-HOM experiment A, L = 5 km.
%
% This diagnostic script solves only the L = 5 km case and plots the mesh,
% velocity components, and pressure in one compact figure.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.5*pi/180;
slope = tan(alpha);
H = 1000;                   % m
bedAmplitude = 500;         % m
L = 5000;                   % m

officialFSMean = 14.20;
officialFSMax = 14.56;

Nx = 10;
Ny = 10;
Nz = 2;
% Nx = 5;
% Ny = 5;
% Nz = 1;

[refnode,elem] = cubemesh([0,1,0,1,0,1],[1/Nx,1/Ny,1/Nz]);
bdFlag = setboundary3(refnode,elem,'Neumann','z==1',...
    'Dirichlet','z==0');

x = L*refnode(:,1);
y = L*refnode(:,2);
zeta = refnode(:,3);
surface = -slope*x;
bed = surface-H+bedAmplitude*sin(2*pi*x/L).*sin(2*pi*y/L);
node = [x,y,bed+zeta.*(surface-bed)];

pde = struct;
pde.A = 1e-16;              % yr^-1 Pa^-3
pde.n = 3;
pde.m = 1;
pde.beta = 0;
pde.rho = rho;
pde.gravity = [0,0,-gravity];
pde.g_N = [];

option.periodic = true;
option.periodic_x = [0,L];
option.periodic_y = [0,L];
option.periodic_slope = [slope,0];
option.pressure_constraint = 'none';
option.eps_reg = 1e-10;
option.maxIt = 150;
option.tol = 1e-8;
option.residual_tol = 1e-8;
option.damping = 0.7;
option.printlevel = 1;
option.quadorder = 4;
option.facequadorder = 4;

[soln,~,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);

[~,edge] = dof3P2(elem);
uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
tolGeometry = 1000*eps(max(1,max(abs(node(:)))));
top = abs(uNode(:,3)+slope*uNode(:,1)) < tolGeometry;
section = top & abs(uNode(:,2)-L/4) < tolGeometry;
if ~any(section)
    topIndex = find(top);
    distance = abs(uNode(topIndex,2)-L/4);
    section = false(size(top));
    section(topIndex(distance == min(distance))) = true;
end

u = soln.ux;
v = soln.uy;
w = soln.uz;
p = soln.p;
speedSection = hypot(u(section),v(section));

fprintf(['ISMIP-HOM A L=5 km: converged=%d, Picard=%3d, ',...
    'mean speed=%9.4e, max speed=%9.4e m/yr\n'],...
    info.converged,info.itStep,mean(speedSection),max(speedSection));
fprintf('Official FS mean=%9.4e, max=%9.4e m/yr\n',...
    officialFSMean,officialFSMax);

figure(1);
set(gcf,'Visible','on');
clf;
plotfields(node,elem,uNode,u,v,w,p,L);
sgtitle('ISMIP-HOM experiment A, L = 5 km');

function plotfields(node,elem,uNode,u,v,w,p,L)
mesh = triangulation(elem,node);
face = freeBoundary(mesh);
nodeValue.u = nodalfield(node,uNode,u);
nodeValue.v = nodalfield(node,uNode,v);
nodeValue.w = nodalfield(node,uNode,w);

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
trisurf(face,node(:,1),node(:,2),node(:,3),...
    'FaceColor','none','EdgeColor',[0.35,0.35,0.35]);
formatfieldaxis(L);
title('mesh');

nexttile;
plotsurfacefield(face,node,nodeValue.u,'velocity u');
formatfieldaxis(L);

nexttile;
plotsurfacefield(face,node,nodeValue.v,'velocity v');
formatfieldaxis(L);

nexttile;
plotsurfacefield(face,node,nodeValue.w,'velocity w');
formatfieldaxis(L);

nexttile;
trisurf(face,node(:,1),node(:,2),node(:,3),p,...
    'FaceColor','interp','EdgeColor','interp');
formatfieldaxis(L);
colorbar;
title('pressure p');

nexttile;
axis off;
end

function nodeValue = nodalfield(node,uNode,value)
nNode = size(node,1);
nodeValue = value(1:nNode);
if size(uNode,1) >= nNode
    tolerance = 1000*eps(max(1,max(abs(node(:)))));
    [tf,loc] = ismembertol(node,uNode,tolerance,'ByRows',true);
    nodeValue(tf) = value(loc(tf));
end
end

function plotsurfacefield(face,node,value,titleText)
trisurf(face,node(:,1),node(:,2),node(:,3),value,...
    'FaceColor','interp','EdgeColor','interp');
colorbar;
title(titleText);
end

function formatfieldaxis(L)
axis equal;
axis tight;
xlim([0,L]);
ylim([0,L]);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(35,22);
end
