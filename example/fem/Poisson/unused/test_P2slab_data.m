function pde = test_P2slab_data
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

    pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, 'g_N',@g_N);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs =  - 2*slope^2 - 2;
    end
    
    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        y1 = y + slope * x;
        u = y1.^2;
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

    % Derivative of the exact solution
    function uprime =  Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = 2*slope*(y + slope*x);
        uprime(:,2) = 2*y + 2*slope*x;
    end

    % g_N = dudn + beta u
    function uN = g_N(p)
        du = Du(p);
        u = exactu(p);
        n = [-slope, -1];
        n = n / norm(n);
        uN = (n(1)*du(1) + n(2)*du(2)) - 1*u;
    end
    
    
end
