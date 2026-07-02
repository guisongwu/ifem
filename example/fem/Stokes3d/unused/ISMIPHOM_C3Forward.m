%% ISMIPHOM_C3FORWARD 3-D ISMIP-HOM experiment C forward solve.
%
% Experiment C is a periodic sloping slab with uniform thickness and
% spatially varying basal friction beta^2(x,y).  The top is traction-free.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.1*pi/180;
slope = tan(alpha);
H = 1000;                   % m

lengthList = 1000*5;        % use 1000*[5;10;20;40;80;160] for a sweep
officialFSMax = 16.00;      % Pattyn et al. 2008, Table 4, L=5 km
officialFSMean = 15.99;     % Pattyn et al. 2008, Table 5, L=5 km

Nx = 4;
Ny = 4;
Nz = 2;

result = table('Size',[numel(lengthList),7],...
    'VariableTypes',{'double','logical','double','double','double',...
                     'double','double'},...
    'VariableNames',{'lengthKm','converged','PicardSteps',...
                     'meanSurfaceUx','maxSurfaceUx',...
                     'officialFSMean','officialFSMax'});

for iLength = 1:numel(lengthList)
    L = lengthList(iLength);
    [node,elem] = cubemesh([0,1,0,1,0,1],[1/Nx,1/Ny,1/Nz]);
    bdFlag = setboundary3(node,elem,'Neumann','z==1',...
        'Robin','z==0');

    % Map the rectangular computational box to the experiment-C sloping
    % slab: z_s=-x*tan(alpha), z_b=z_s-H.
    node(:,1) = L*node(:,1);
    node(:,2) = L*node(:,2);
    zeta = node(:,3);
    node(:,3) = -node(:,1)*slope-H+zeta*H;

    pde = struct;
    pde.A = 1e-16;              % yr^-1 Pa^-3
    pde.n = 3;
    pde.m = 1;                  % linear sliding
    pde.beta = @(pt) 1000+1000*sin(2*pi*pt(:,1)/L).*...
        sin(2*pi*pt(:,2)/L);
    pde.rho = rho;
    pde.gravity = [0,0,-gravity];
    pde.g_N = [];

    option.periodic = true;
    option.periodic_x = [0,L];
    option.periodic_y = [0,L];
    option.periodic_slope_x = slope;
    option.pressure_constraint = 'none';
    option.eps_reg = 1e-10;
    option.maxIt = 120;
    option.tol = 1e-8;
    option.residual_tol = 1e-8;
    option.damping = 0.7;
    option.printlevel = 1;
    option.quadorder = 4;
    option.facequadorder = 4;

    [soln,eqn,info] = NonlinearStokes3P2P1(node,elem,bdFlag,pde,option);
    [~,edge] = dof3P2(elem);
    uNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
    top = find(abs(uNode(:,3)+slope*uNode(:,1)) < ...
        1000*eps(max(1,max(abs(node(:))))));

    result.lengthKm(iLength) = L/1000;
    result.converged(iLength) = info.converged;
    result.PicardSteps(iLength) = info.itStep;
    result.meanSurfaceUx(iLength) = mean(soln.ux(top));
    result.maxSurfaceUx(iLength) = max(soln.ux(top));
    result.officialFSMean(iLength) = officialFSMean(min(iLength,end));
    result.officialFSMax(iLength) = officialFSMax(min(iLength,end));

    fprintf(['ISMIP-HOM C3 L=%6.1f km: converged=%d, Picard=%3d, ',...
        'mean ux=%9.4e, max ux=%9.4e m/yr\n'],...
        L/1000,info.converged,info.itStep,...
        result.meanSurfaceUx(iLength),result.maxSurfaceUx(iLength));
end

disp(result);
