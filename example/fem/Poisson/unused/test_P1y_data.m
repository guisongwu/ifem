function pde = test_P1y_data
%% SINCOSDATA trigonometric  data for Poisson equation
%
%     f = 0;
%     u = x - y;
%     Du = (1, -1);
%
% TODO: Neumann and Robin
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

    pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs =  0;
    end
    
    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        u = y;
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

    % Derivative of the exact solution
    function uprime =  Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = 0;
        uprime(:,2) = 1;
    end

    
end
