%% ISMIPHOM_B Flowline version of the ISMIP-HOM experiment B.
%
% Experiment B is the two-dimensional flowline version of ISMIP-HOM A:
% flow over a sinusoidal bed with no slip at the base.  The surface is a
% plane slope and the left/right sides are periodic.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.5*pi/180;         % surface slope for ISMIP-HOM A/B
slope = tan(alpha);
meanThickness = 1000;       % m
bedAmplitude = 500;         % m

% Standard ISMIP-HOM length scales.  Use one entry here for a faster smoke
% test, or all entries to reproduce the diagnostic length-scale sweep.
lengthList = 1000*[5;10;20;40;80;160];
officialFSMax = [11.76;22.82;46.91;73.77;95.12;108.33];
officialFSMean = [11.04;19.09;28.28;35.75;39.76;41.40];
officialFSCurveData = readofficialfscurvedata(...
    'HOM_B_official_fs_curve.csv');

% Keep the same number of cells for each normalized domain.  The physical
% mesh spacing therefore grows with L, as in ISMIPHOM_D.m.
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
    [xSurface,uxSurface] = getsurfacevelocity(...
        node,elem,soln,L,slope);

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

    fprintf(['ISMIP-HOM B L=%6.1f km: converged=%d, Picard=%3d, ',...
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
    yline(officialFSMean(iLength),'b--','LineWidth',1.1,...
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
sgtitle('ISMIP-HOM experiment B: surface horizontal velocity');

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

function [xSurface,uxSurface] = getsurfacevelocity(node,elem,soln,L,slope)
    [~,edge] = dofP2(elem);
    velocityNode = [node;(node(edge(:,1),:)+node(edge(:,2),:))/2];
    xWrapped = mod(velocityNode(:,1),L);
    surface = -slope*velocityNode(:,1);
    tolerance = 100*eps(max(1,max(abs(velocityNode(:)))));
    topDof = find(abs(velocityNode(:,2)-surface) <= tolerance);
    xSurface = xWrapped(topDof);
    uxSurface = soln.ux(topDof);
    [xSurface,order] = sort(xSurface);
    uxSurface = uxSurface(order);
end
