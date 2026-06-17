function pde = stokes_data_P4
%% STOKESDATA2 data for Stokes equations
%
% The solution u is a polynomial satisfying the zero Dirichlet boundary
% condition. 
%
% Created by Ming Wang.
%
    pde = struct('f', @f, 'exactp', @exactp, ...
                 'exactu',@exactu, ...
                 'exactux',@exactux, ...
                 'exactuy',@exactuy, ...
                 'g_D', @exactu, 'g_N', @g_N, 'g_R',@g_R, 'g_RN', @g_RN, 'g_Dn', @g_Dn);
%

global slope;

%% subfunction

% x*power(y,2) + power(y,3),
% power(x,3) - power(y,3)/3.,
% x*y - power(y,2),
% -2*x - 5*y,
% -5*x

    function z = exactux(pt)
        x = pt(:,1); y = pt(:,2);
        z = x.*power(y,2) + power(y,3);
    end
    function z = exactuy(pt)
        x = pt(:,1); y = pt(:,2);
        z = power(x,3) - power(y,3)/3.;
    end


    function z = exactu(pt)
        x = pt(:,1); y = pt(:,2);
        z(:,1) = exactux(pt);
        z(:,2) = exactuy(pt);
    end

    function z = exactp(pt)
        x = pt(:,1); y = pt(:,2);
        z = x.*y - power(y,2);
    end

    function z = f(pt) % load data (right hand side function)
        x = pt(:,1); y = pt(:,2);
        z(:,1) = x*0 -2*x - 5*y;
        z(:,2) = x*0 -5*x;
    end


    % power(y,2)
    % 2*x*y + 3*power(y,2)
    % 3*power(x,2)
    % -power(y,2)
    % y
    % x - 2*y
    
    % Derivative of the exact solution
    function du =  Du(pt)
        x = pt(:,1); y = pt(:,2);
        du(:,1,1) = x*0 + power(y,2);
        du(:,1,2) = x*0 + 2*x.*y + 3*power(y,2);
        du(:,2,1) = x*0 + 3*power(x,2);
        du(:,2,2) = x*0 -power(y,2);
    end

    function dp =  Dp(pt)
        x = pt(:,1); y = pt(:,2);
        dp(:,1) = x*0 + y;
        dp(:,2) = x*0 + x - 2*y;
    end




    % g_N = dudn
    function uN = g_N(pt)
        x = pt(:,1); y = pt(:,2);
        y1 = y + slope * x;
        du = Du(pt);
        p = exactp(pt);

        n = [slope, 1];
        n = n / norm(n);

        % size(du)
        % size(p)
        
        uN(:,1) = n(1)*du(:,1,1) + n(2)*du(:,1,2) - n(1) * p(:);
        uN(:,2) = n(1)*du(:,2,1) + n(2)*du(:,2,2) - n(2) * p(:);
    end



    % g_RN = dudn + beta u
    function uR = g_R(pt)
        x = pt(:,1); y = pt(:,2);
        uR = 0*x + 1.+0.1*x;
    end
    
    function uN = g_RN(pt)
        x = pt(:,1); y = pt(:,2);
        y1 = y + slope * x;
        du = Du(pt);

        u = exactu(pt);
        p = exactp(pt);
        du = Du(pt);

        n = [-slope, -1];
        n = n / norm(n);

        beta = g_R(pt);

        uN(:,1) = n(1)*du(:,1,1) + n(2)*du(:,1,2) - n(1) * p(:) + beta .* u(:, 1);
        uN(:,2) = n(1)*du(:,2,1) + n(2)*du(:,2,2) - n(2) * p(:) + beta .* u(:, 2);
    end


    function uN = g_Dn(pt)
        x = pt(:,1); y = pt(:,2);
        y1 = y + slope * x;

        u = exactu(pt);
        n = [-slope, -1];
        n = n / norm(n);

        uN = u(:,1) * n(1) + u(:,2) * n(2);
    end

end
