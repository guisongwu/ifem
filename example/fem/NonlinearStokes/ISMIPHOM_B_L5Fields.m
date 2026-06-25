%% ISMIPHOM_B_L5FIELDS Plot fields for ISMIP-HOM experiment B at L = 5 km.
%
% This diagnostic script solves only the L = 5 km flowline case and plots
% the mesh, horizontal velocity, vertical velocity, and pressure in one
% 2-by-2 figure.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.5*pi/180;         % surface slope for ISMIP-HOM A/B
slope = tan(alpha);
meanThickness = 1000;       % m
bedAmplitude = 500;         % m
L = 5000;                   % m

Nx = 40;
Nz = 10;

[node,elem] = rectanglemesh(L,1,Nx,Nz);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node = maptoexperimentb(node,L,slope,meanThickness,bedAmplitude);

pde = struct;
pde.A = 1e-16;              % yr^-1 Pa^-3
pde.n = 3;
pde.m = 1;
pde.beta = 0;               % unused for option.bed_condition = 'no-slip'
pde.rho = rho;
pde.gravity = [0,-gravity];
pde.g_N = [];               % traction-free upper surface

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-10;
option.maxIt = 150;
option.tol = 1e-8;
option.residual_tol = 1e-8;
option.damping = 0.7;
option.bed_condition = 'no-slip';
option.printlevel = 0;
option.quadorder = 4;

[soln,~,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);
fprintf('ISMIP-HOM B L=5 km: converged=%d, Picard=%d\n',...
    info.converged,info.itStep);

figure(1);
set(gcf,'Visible','on');
clf;
plotfields(node,elem,soln);
sgtitle('ISMIP-HOM experiment B, L = 5 km');

function node = maptoexperimentb(node,L,slope,meanThickness,bedAmplitude)
    x = node(:,1);
    q = node(:,2);
    surface = -slope*x;
    bed = surface-meanThickness+bedAmplitude*sin(2*pi*x/L);
    node(:,2) = bed+q.*(surface-bed);
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
            n1 = cellId(iz,ix);
            n2 = cellId(iz,ix+1);
            n3 = cellId(iz+1,ix+1);
            n4 = cellId(iz+1,ix);
            cursor = cursor+1;
            elem(cursor,:) = [n2,n3,n1];
            cursor = cursor+1;
            elem(cursor,:) = [n4,n1,n3];
        end
    end
end

function plotfields(node,elem,soln)
    nNode = size(node,1);

    subplot(2,2,1);
    showmesh(node,elem);
    axis equal;
    axis tight;
    title('mesh');

    subplot(2,2,2);
    trisurf(elem,node(:,1),node(:,2),soln.ux(1:nNode),...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title('horizontal velocity u_x');
    view(2);

    subplot(2,2,3);
    trisurf(elem,node(:,1),node(:,2),soln.uz(1:nNode),...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title('vertical velocity u_z');
    view(2);

    subplot(2,2,4);
    trisurf(elem,node(:,1),node(:,2),soln.p,...
        'FaceColor','interp','EdgeColor','interp');
    axis equal;
    axis tight;
    colorbar;
    title('pressure p');
    view(2);
end
