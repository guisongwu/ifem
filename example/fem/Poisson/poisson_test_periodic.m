
%% Test Poisson

close all; 
clear variables;

global slope;                           % slab slope
slope = 0.1;

option.verb = 1;
option.gNquadorder = 4;
option.periodic = false;


figure(1);
%set(gcf, 'Position', get(0, 'Screensize'));
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);    

maxIt = 3; 
errL2 = zeros(maxIt,1); 



h = 0.1;

for k = 1:maxIt


    % -------------
    % Set test case
    % -------------
    %pde = test_P1_data;
    %pde = test_P1slab_data;
    %pde = test_P2slab_data;
    pde = poisson_data_P2_slab;
    %pde = test_P2_data;


    
    [node,elem] = squaremesh([0 1 0 1], h);
    bdFlag = setboundary(node,elem, 'Neumann','y==1', 'Robin', 'y==0');
    %bdFlag = setboundary(node,elem, 'Dirichlet','y==0|y==1');

    
    
    if slope ~= 0
        %  
        % Slab test remap coord
        % Note: this must be done after bdFlag is set
        % 
        fprintf(2, 'Slab test %f\n', slope);
        node(:,2) = node(:,2) - slope * node(:,1);
    end

    if true
        % 
        %  
        % Change func to data
        % 
        %
        warning('\nChange function to data !!!\n');
        
        xbot = [0:h:1-h]';                     % periodic without last one
        ybot = 0 - slope * xbot;
        pbot = [xbot ybot];

        ytop = 1 - slope * xbot;
        ptop = [xbot ytop];
        
        pde.g_N = pde.g_N(ptop);
        pde.g_R = pde.g_R(pbot);
        pde.g_RN = pde.g_RN(pbot);

        % pde.g_N
        % pde.g_R
        % pde.g_RN

        pde
    else
        pde
    end

    
    [soln,eqn,info] = Poisson_periodic(node,elem,bdFlag,pde,option);
    uh = soln.u;


    if false
        % 
        % Check res
        % 
        uE = pde.exactu(node);
        up = zeros(length(unique(eqn.Ip)), 1);
        for i = 1:length(node)
            up(eqn.Ip(i)) = uE(i);
        end
        fprintf('res: \n');
        eqn.Ap * up - eqn.bp
    end
    
    showresult(node,elem,uh);

    errL2(k) = getL2error(node,elem,pde.exactu,uh);
    fprintf('pass: %d, h: %f, err: %e\n', k, h, errL2(k));

    h = h / 2;
end


if maxIt > 1
    fprintf('\nConvergence rate\n');
    for k = 1:maxIt
        if k == 1
            fprintf('err: %e -\n', errL2(k));
        else
            fprintf('err: %e %.2f\n', errL2(k), -log(errL2(k)/errL2(k-1)) / log(2) );
        end
    end
end
