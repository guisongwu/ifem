function pde = test_P1slab_data
%% SINCOSDATA trigonometric  data for Poisson equation
%
%     f = 0;
%     u = x - y;
%     Du = (1, -1);
%
% TODO: Neumann and Robin
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.
    global slope;

    pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, 'g_N',@g_N, 'g_R',@g_R, 'g_RN',@g_RN);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs = 0*x + 0;
    end
    
    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        y1 = y + slope * x;
        u = y1;
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

    % Derivative of the exact solution
    function uprime =  Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = x*0 + slope;
        uprime(:,2) = x*0 + 1;
    end

    function uR = g_R(p)
        x = p(:,1); y = p(:,2);
        uR = 0*x + 1.;
    end

    
    % g_N = dudn
    function uN = g_N(p)
        x = p(:,1); y = p(:,2);
        y1 = y + slope * x;
        du = Du(p);

        n = [slope, 1];
        n = n / norm(n);
        uN = n(1)*du(:,1) + n(2)*du(:,2);
    end

    % g_RN = dudn + beta u
    function uN = g_RN(p)
        x = p(:,1); y = p(:,2);
        y1 = y + slope * x;
        du = Du(p);

        u = exactu(p);
        n = [-slope, -1];
        n = n / norm(n);
        beta = g_R(p);
        uN = (n(1)*du(:,1) + n(2)*du(:,2)) + beta.*u;
    end
    
    
end
