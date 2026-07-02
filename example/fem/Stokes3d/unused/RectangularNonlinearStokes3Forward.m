%% RECTANGULARNONLINEARSTOKES3FORWARD Nonlinear Stokes solve in a cuboid.
%
% This is a 3-D counterpart of the NonlinearStokes forward examples.  The
% side walls are no-slip, the top is traction-free, and the flat bed uses a
% nonlinear Weertman sliding law.

close all;
clear variables;

%% Geometry and mesh
L = 2;
W = 1;
H = 1;
h = 0.5;
[node,elem] = cubemesh([0,L,0,W,0,H],h);
bdFlag = setboundary3(node,elem,'Dirichlet','all',...
    'Neumann',sprintf('z==%.17g',H),'Robin','z==0');

%% Nonlinear Stokes model
pde = struct;
pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.rho = 1;
pde.gravity = [0.15,0,-1];
pde.beta = @(pt) 2*(1+0.25*cos(2*pi*pt(:,1)/L).*...
    cos(2*pi*pt(:,2)/W));
pde.g_N = @(pt) zeros(size(pt,1),3);

option.eps_reg = 1e-3;
option.maxIt = 80;
option.tol = 1e-9;
option.residual_tol = 1e-9;
option.damping = 0.8;
option.printlevel = 1;
option.quadorder = 4;
option.facequadorder = 4;
option.assemble_tangent = false;

[soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);
fprintf('Converged: %d, Picard steps: %d, final residual: %.04e\n',...
    info.converged,info.itStep,info.nonlinearResidual(end));

%% Simple visualization on vertex values
N = size(node,1);
speed = sqrt(soln.ux(1:N).^2+soln.uy(1:N).^2+soln.uz(1:N).^2);

figure(1);
scatter3(node(:,1),node(:,2),node(:,3),40,speed,'filled');
axis equal tight;
view(3);
colorbar;
title('P2 velocity speed at vertices');

figure(2);
scatter3(node(:,1),node(:,2),node(:,3),40,soln.p,'filled');
axis equal tight;
view(3);
colorbar;
title('P1 pressure');
