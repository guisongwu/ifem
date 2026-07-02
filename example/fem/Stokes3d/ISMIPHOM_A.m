%% ISMIPHOM_A 3-D ISMIP-HOM experiment A forward solve.
%
% Experiment A uses a periodic sloping slab with sinusoidal bedrock and a
% stress-free upper surface.  The basal boundary is no slip.

close all;
clear variables;

rho = 910;                  % kg/m^3
gravity = 9.81;             % m/s^2
alpha = 0.5*pi/180;
slope = tan(alpha);
H = 1000;                   % m
bedAmplitude = 500;         % m

lengthList = 1000*[5;10;20;40;80;160];
officialLength = [5;10;20;40;80;160];
officialFSMax = [14.56;24.36;39.73;63.89;87.10;102.63];
officialFSMean = [14.20;20.02;24.74;31.89;37.31;39.98];

Nx = 10;
Ny = 10;
Nz = 2;

result = table('Size',[numel(lengthList),7],...
    'VariableTypes',{'double','logical','double','double','double',...
                     'double','double'},...
    'VariableNames',{'lengthKm','converged','PicardSteps',...
                     'meanSectionSpeed','maxSectionSpeed',...
                     'officialFSMean','officialFSMax'});
surfaceProfile = cell(numel(lengthList),1);

for iLength = 1:numel(lengthList)
    L = lengthList(iLength);
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
    speed = hypot(soln.ux(section),soln.uy(section));
    xSection = mod(uNode(section,1),L);
    [xSection,order] = sort(xSection);
    speedProfile = speed(order);

    result.lengthKm(iLength) = L/1000;
    result.converged(iLength) = info.converged;
    result.PicardSteps(iLength) = info.itStep;
    result.meanSectionSpeed(iLength) = mean(speed);
    result.maxSectionSpeed(iLength) = max(speed);
    officialIndex = find(officialLength == L/1000,1);
    if ~isempty(officialIndex)
        result.officialFSMean(iLength) = officialFSMean(officialIndex);
        result.officialFSMax(iLength) = officialFSMax(officialIndex);
    else
        result.officialFSMean(iLength) = NaN;
        result.officialFSMax(iLength) = NaN;
    end
    surfaceProfile{iLength} = [xSection/L,speedProfile];

    fprintf(['ISMIP-HOM A L=%6.1f km: converged=%d, Picard=%3d, ',...
        'mean speed=%9.4e, max speed=%9.4e m/yr\n'],...
        L/1000,info.converged,info.itStep,...
        result.meanSectionSpeed(iLength),result.maxSectionSpeed(iLength));
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
    yline(officialFSMax(iLength),'k--','LineWidth',1.1,...
        'DisplayName','official FS max');
    yline(officialFSMean(iLength),'b:','LineWidth',1.1,...
        'DisplayName','official FS mean');
    hold off;
    grid on;
    xlabel('x/L');
    ylabel('surface speed (m/yr)');
    title(sprintf('L = %.0f km',lengthList(iLength)/1000));
    if iLength == 1
        legend('Location','best');
    end
end
sgtitle('ISMIP-HOM experiment A: surface speed at y = L/4');
