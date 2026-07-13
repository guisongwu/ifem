%% ISMIPHOM_D Flowline version of the ISMIP-HOM experiment D.
%
% Experiment D is the two-dimensional flowline version of ISMIP-HOM C:
% flow over a flat sloping bed with a spatially varying linear basal
% friction coefficient beta^2.  The problem is diagnostic: geometry is
% fixed and only the Stokes velocity and pressure are solved.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.1*pi/180;         % surface and bed slope for ISMIP-HOM C/D
slope = tan(alpha);
H = 1000;                   % m

% Standard ISMIP-HOM length scales.  Use one entry here for a faster smoke
% test, or all entries to reproduce the diagnostic length-scale sweep.
lengthList = 1000*[5;10;20;40;80;160];
officialFSMax = [16.48;17.11;21.33;41.51;97.64;238.44];
officialFSMean = [16.43;16.81;18.40;24.63;37.00;57.17];
officialFSCurveData = readofficialfscurvedata(...
    'HOM_D_official_fs_curve.csv');

% Keep the vertical resolution fixed and scale only the horizontal spacing.
% This avoids excessively large meshes for the 160 km case.
Nx = 40;
Nz = 10;

result = table('Size',[numel(lengthList),8],...
    'VariableTypes',{'double','logical','double','double','double',...
                     'double','double','double'},...
    'VariableNames',{'lengthKm','converged','PicardSteps',...
                     'meanSurfaceUx','maxSurfaceUx',...
                     'officialFSMean','officialFSMax','maxDifference'});
surfaceProfile = cell(numel(lengthList),1);
officialFSCurveProfile = cell(numel(lengthList),1);

for iLength = 1:numel(lengthList)
    L = lengthList(iLength);
    [node,elem] = rectanglemesh(L,H,Nx,Nz);

    topBoundaryExpression = sprintf('y==%.17g',H);
    bdFlag = setboundary(node,elem,'Neumann',topBoundaryExpression,...
        'Robin','y==0');
    node(:,2) = node(:,2)-slope*node(:,1);

    pde = struct;
    pde.A = 1e-16;              % yr^-1 Pa^-3
    pde.n = 3;
    pde.m = 1;                  % ISMIP-HOM uses linear sliding
    pde.beta = @(pt) ismiphombeta(pt,L);
    pde.rho = rho;
    pde.gravity = [0,-gravity];  % body force rho*g has units Pa/m
    pde.g_N = [];               % stress-free upper surface

    option.periodic = true;
    option.periodic_x = [0,L];
    option.eps_reg = 1e-10;
    option.maxIt = 150;
    option.tol = 1e-8;
    option.residual_tol = 1e-8;
    option.damping = 0.7;
    option.printlevel = 0;
    option.quadorder = 4;

    [soln,~,info] = NonlinearStokesP2P1(node,elem,bdFlag,pde,option);
    [xSurface,uxSurface] = getsurfacevelocity(node,elem,soln,L,H,slope);

    result.lengthKm(iLength) = L/1000;
    result.converged(iLength) = info.converged;
    result.PicardSteps(iLength) = info.itStep;
    result.meanSurfaceUx(iLength) = mean(uxSurface);
    result.maxSurfaceUx(iLength) = max(uxSurface);
    result.officialFSMean(iLength) = officialFSMean(iLength);
    result.officialFSMax(iLength) = officialFSMax(iLength);
    officialFSCurveProfile{iLength} = officialFSCurveData{iLength};
    result.maxDifference(iLength) = max(uxSurface)-officialFSMax(iLength);
    surfaceProfile{iLength} = [xSurface/L,uxSurface];

    fprintf(['ISMIP-HOM D L=%6.1f km: converged=%d, Picard=%3d, ',...
        'max surface ux=%9.4e m/yr\n'],...
        L/1000,info.converged,info.itStep,max(uxSurface));
end

disp(result);

figure(1);
set(gcf,'Visible','on');
clf;
for iLength = 1:numel(lengthList)
    subplot(2,3,iLength);
    profile = surfaceProfile{iLength};
    plot(profile(:,1),profile(:,2),'r-','LineWidth',1.6,...
        'DisplayName','this solver');
    hold on;
    if ~isempty(officialFSCurveProfile{iLength})
        fsProfile = officialFSCurveProfile{iLength};
        plot(fsProfile(:,1),fsProfile(:,2),'-.','LineWidth',1.5,...
            'Color',[0.85,0.33,0.10],...
            'DisplayName','official FS mean curve');
    end
    yline(officialFSMax(iLength),'k--','LineWidth',1.1,...
        'DisplayName','official FS max');
    yline(officialFSMean(iLength),'b:','LineWidth',1.1,...
        'DisplayName','official FS mean');
    hold off;
    grid on;
    xlabel('x/L');
    ylabel('surface u_x (m/yr)');
    title(sprintf('L = %.0f km',lengthList(iLength)/1000));
    if iLength == 1
        legend('Location','best');
    end
end
sgtitle('ISMIP-HOM experiment D: surface horizontal velocity');

function beta = ismiphombeta(pt,L)
    x = mod(pt(:,1),L);
    % Pattyn et al. define the basal friction coefficient beta^2 directly:
    % beta^2(x) = 1000+1000*sin(2*pi*x/L), in Pa yr m^{-1}.
    % NonlinearStokesP2P1 uses pde.beta as the coefficient multiplying
    % tangential velocity, so do not square this value again.
    beta = 1000+1000*sin(2*pi*x/L);
end

function officialFSCurveData = readofficialfscurvedata(filename)
    lengthKm = [5;10;20;40;80;160];
    officialFSCurveData = cell(numel(lengthKm),1);
    if ~isfile(filename)
        return;
    end

    data = readtable(filename);
    for iLength = 1:numel(lengthKm)
        isCurrentLength = data.lengthKm == lengthKm(iLength);
        officialFSCurveData{iLength} = [data.xOverL(isCurrentLength),...
            data.uxFSMean(isCurrentLength)];
    end
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

function [xSurface,uxSurface] = getsurfacevelocity(node,elem,soln,L,H,slope)
    [~,edge] = dofP2(elem);
    velocityNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
    q = velocityNode(:,2)+slope*velocityNode(:,1);
    tolerance = 100*eps(max(1,max(abs(velocityNode(:)))));
    topDof = find(abs(q-H) <= tolerance);
    xSurface = mod(velocityNode(topDof,1),L);
    uxSurface = soln.ux(topDof);
    [xSurface,order] = sort(xSurface);
    uxSurface = uxSurface(order);
end
