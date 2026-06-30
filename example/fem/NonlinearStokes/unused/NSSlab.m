%% NSSLAB Nonlinear full-Stokes flow in an ice slab.
%
% The physical coordinates are (x,z).  The domain is
%
%   0 < x < L,  -s*x < z < H-s*x.
%
% Left and right boundaries are periodic, the upper surface is
% traction-free, and the bed is impermeable with a Weertman sliding law.

close all;
clear variables;

% L = 1;
% H = 0.5;
L = 5;
H = 1;

slope = 0.1;
h = 1/10;
% betaCase = 'constant';
% betaCase = 'linear';
% betaCase = 'quadratic';
betaCase = 'sin';
% betaCase = 'cos';
% betaCase = 'bump';

[node,elem] = squaremesh([0,L,0,H],h);
topBoundaryExpression = sprintf('y==%.17g',H);
bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
    'Robin','y==0');
node(:,2) = node(:,2)-slope*node(:,1);

% Dimensionless demonstration parameters.  Replace these by consistently
% scaled physical values when working in SI or glaciological units.
pde = struct;
pde.A = 1;
pde.n = 3;
pde.beta = makebedbeta(betaCase,L);
pde.m = 1/3;
pde.rho = 1;
pde.gravity = [0,-1];
pde.g_N = [];                 % traction-free upper surface

option.periodic = true;
option.periodic_x = [0,L];
option.eps_reg = 1e-9;
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
set(gcf,'Visible','on');
clf;
subplot(2,2,1);
showmesh(node,elem);
axis equal;
axis tight;
title('mesh');

subplot(2,2,2);
trisurf(elem,node(:,1),node(:,2),soln.ux(1:size(node,1)),...
    'FaceColor','interp','EdgeColor','interp');
axis equal;
axis tight;
colorbar;
title('horizontal velocity u_x');
view(2);

subplot(2,2,3);
trisurf(elem,node(:,1),node(:,2),soln.uz(1:size(node,1)),...
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
title('pressure');
view(2);

figure(2);
xPlot = linspace(0,L,401)';
bedPoint = [xPlot,-slope*xPlot];
betaPlot = evaluatebeta(pde.beta,bedPoint);
plot(xPlot,betaPlot,'k-','LineWidth',1.6);
grid on;
xlabel('x');
ylabel('\beta');
title(sprintf('bed beta: %s',betaCase));

function beta = makebedbeta(betaCase,L)
    switch lower(betaCase)
        case 'constant'
            beta = 10;
        case 'linear'
            beta = @(pt) 8+4*periodicx(pt(:,1),L)/L;
        case 'quadratic'
            x = @(pt) periodicx(pt(:,1),L)/L;
            beta = @(pt) 8+8*(x(pt)-0.5).^2;
        case 'sin'
            beta = @(pt) 10+2*sin(2*pi*periodicx(pt(:,1),L)/L);
        case 'cos'
            beta = @(pt) 10+2*cos(2*pi*periodicx(pt(:,1),L)/L);
        case 'bump'
            x = @(pt) periodicx(pt(:,1),L)/L;
            beta = @(pt) 8+5*exp(-((x(pt)-0.5)/0.18).^2);
        otherwise
            error('Unknown betaCase: %s.',betaCase);
    end
end

function xWrapped = periodicx(x,L)
    xWrapped = mod(x,L);
end

function value = evaluatebeta(beta,pt)
    if isa(beta,'function_handle')
        value = beta(pt);
    else
        value = beta+zeros(size(pt,1),1);
    end
end
