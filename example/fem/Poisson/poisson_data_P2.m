function pde = test_P2_data
%% SINCOSDATA trigonometric  data for Poisson equation
%
%     f = -4;
%     u = (x + y)^2
%     Du = (x+y, x+y);
%
% TODO: Neumann and Robin
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

    pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs =  -4;
    end
    
    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        u = (x + y).^2;
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

    % Derivative of the exact solution
    function uprime =  Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = x+y;
        uprime(:,2) = x+y;
    end

    
end
