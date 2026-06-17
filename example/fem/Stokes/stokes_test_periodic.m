
%% Test Stokes periodic

close all; 
clear variables;

global slope dbg_case h;                           % slab slope
slope = 0.1;

option.verb = 1;
option.gNquadorder = 4;
option.solver = 'direct';
option.use_slip = true;                 % robin -> slip

test_periodic = true;

figure(1);
%set(gcf, 'Position', get(0, 'Screensize'));
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);    

maxIt = 3; 
errL2 = zeros(maxIt,1); 

dbg_case = 0;

h = 0.1;

for k = 1:maxIt


    % -------------
    % Set test case
    % -------------
    if ~test_periodic
        %pde = stokes_data_P4;
        %pde = stokes_data_P1;
        pde = stokes_data_P2_slab;
        %pde = stokes_data_P1;
        %pde = stokes_data_P2_period;
        option.periodic = false;
    else
        %pde = stokes_data_P2_period;
        pde = stokes_data_sin_period;
        option.periodic = true;
        warning('Test periodic');
    end

    
    [node,elem] = squaremesh([0 1 0 1], h);
    %bdFlag = setboundary(node,elem, 'Dirichlet');
    %bdFlag = setboundary(node,elem, 'Dirichlet','~(y==1)', 'Neumann', 'y==1');
    %bdFlag = setboundary(node,elem, 'Dirichlet','x==0 | x==1', 'Neumann', 'y==1', 'Robin', 'y==0');
    bdFlag = setboundary(node,elem, 'Neumann','y==1', 'Robin', 'y==0');

    
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
        
        xbot = [0:h:1-h h/2:h:1]';                     % periodic without last one
        ybot = 0 - slope * xbot;
        pt_bot = [xbot ybot];

        ytop = 1 - slope * xbot;
        pt_top = [xbot ytop];
        
        pde.g_N = pde.g_N(pt_top);
        pde.g_R = linearize_bot(pde.g_R(pt_bot));
        %pde.g_R = pde.g_R(pt_bot);
        pde.g_RN = pde.g_RN(pt_bot);
        pde.g_Dn = pde.g_Dn(pt_bot);

        % pde.g_N
        % pde.g_R
        % pde.g_RN
    else
    end

    
    [soln,eqn,info] = StokesP2P1_periodic(node,elem,bdFlag,pde,option);
    uh = soln.u;
    ph = soln.p;

    Nu = length(uh(:)) / 2;

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
    
    if length(uh)+length(ph) < 2e3
        figure(1);  showresult(node,elem,ph);    
        figure(2);  showresult(node,elem,uh(1:Nu));    
        figure(3);  showresult(node,elem,uh(Nu+1:2*Nu));    
    end


    % error
    uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);

    % error u H1 norm using A
    % erru(k) = sqrt((uh-uI(:))'*eqn.A*(uh-uI(:)));

    % error u L2 norm
    eux = getL2error(node,elem,pde.exactux,uh(1:Nu));
    euy = getL2error(node,elem,pde.exactuy,uh(Nu+1:2*Nu));
    erru(k) = norm([eux, euy]);

    errp(k) = getL2error(node,elem,pde.exactp,ph);
    fprintf('pass: %d, h: %f, err: %e %e\n', k, h, erru(k), errp(k));
    
    h = h / 2;
end


if maxIt > 1
    fprintf('\nConvergence rate u\n');
    for k = 1:maxIt
        if k == 1
            fprintf('erru: %e -\n', erru(k));
        else
            fprintf('erru: %e %.2f\n', erru(k), -log(erru(k)/erru(k-1)) / log(2) );
        end
    end
    fprintf('\nConvergence rate p\n');
    for k = 1:maxIt
        if k == 1
            fprintf('errp: %e -\n', errp(k));
        else
            fprintf('errp: %e %.2f\n', errp(k), -log(errp(k)/errp(k-1)) / log(2) );
        end
    end
end
