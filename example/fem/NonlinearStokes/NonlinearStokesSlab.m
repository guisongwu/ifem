%% NONLINEARSTOKESSLAB Nonlinear full-Stokes flow in an ice slab.
%
% The physical coordinates are (x,z).  The domain is
%
%   0 < x < L,  -s*x < z < H-s*x.
%
% Left and right boundaries are periodic, the upper surface is
% traction-free, and the bed is impermeable with a Weertman sliding law.

close all;
clear variables;

L = 1;
H = 1;
slope = 0.1;
h = 1/8;

[node,elem] = squaremesh([0,L,0,H],h);
bdFlag = setboundary(node,elem,'Neumann','y==1','Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

% Dimensionless demonstration parameters.  Replace these by consistently
% scaled physical values when working in SI or glaciological units.
pde = struct;
pde.A = 1;
pde.n = 3;
pde.beta = 10;
pde.m = 1/3;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];                 % traction-free upper surface

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-4;
option.maxIt = 100;
option.tol = 1e-8;
option.damping = 0.8;
option.printlevel = 1;

[soln,eqn,info] = NonlinearStokesP2P1(...
    node,elem,bdFlag,pde,option);

Nu = length(soln.ux);
fprintf('Converged: %d, nonlinear iterations: %d\n',...
    info.converged,info.itStep);

figure(1);
showresult(node,elem,soln.p);
title('pressure');

figure(2);
showresult(node,elem,soln.ux(1:size(node,1)));
title('horizontal velocity u_x');

figure(3);
showresult(node,elem,soln.uz(1:size(node,1)));
title('vertical velocity u_z');
