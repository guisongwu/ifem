function pde = test_Stokes_P1_data
%% STOKESDATA2 data for Stokes equations
%
% The solution u is a polynomial satisfying the zero Dirichlet boundary
% condition. 
%
% Created by Ming Wang.
%
    pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu); 
    %
    %% subfunction


    function z = f(p) % load data (right hand side function)
        x = p(:,1); y = p(:,2);
        z(:,1) = 0;
        z(:,2) = 0;
    end

    function z = exactu(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = x + y;
        z(:,2) = x - y;
    end

    function z = exactp(p)
        x = p(:,1); y = p(:,2);
        z = 0; %2*x;
    end

end