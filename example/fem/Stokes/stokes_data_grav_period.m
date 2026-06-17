function pde = stokes_data_grav_period
%% STOKESDATA2 data for Stokes equations
%
% Flow under gravity, when m0 is const, we have exact soltuion
%
% Created by Ming Wang.
%
    pde = struct('f', @f, ...
                 'fp', @fp, ...
                 'exactp', @exactp, ...
                 'exactu',@exactu, ...
                 'exactux',@exactux, ...
                 'exactuy',@exactuy, ...
                 'g_D', @exactu, ...
                 'g_N', @g_N, ...
                 'g_R',@g_R, ...
                 'g_RN', @g_RN, ...
                 'g_Dn', @g_Dn);
%

global slope;


hh = 0.5;
m0 = 1;




%% subfunction
    function z = exactux(pt)
        x = pt(:,1); y = pt(:,2);
        yy = y + x * slope;
        z = 0*x +  1/m0 * hh / (1 + slope^2) - slope/2/(1 + slope^2)^2 * yy .* (yy - 2 * hh); 
    end
    function z = exactuy(pt)
        x = pt(:,1); y = pt(:,2);
        z = -slope * exactux(pt);
    end

    function z = exactu(pt)
        x = pt(:,1); y = pt(:,2);
        z(:,1) = exactux(pt);
        z(:,2) = exactuy(pt);
    end

    function z = exactp(pt)
        x = pt(:,1); y = pt(:,2);
        yy = y + x * slope;
        z = 0*x + 1 / (1 + slope^2) * (hh - yy);
    end

    function z = f(pt) % load data (right hand side function)
        x = pt(:,1); y = pt(:,2);
        z(:,1) = 0*x;
        z(:,2) = 0*x - 1;
    end

    function z = fp(pt) % load data (right hand side function)
        x = pt(:,1); y = pt(:,2);
        z = 0*x;
    end




    % Derivative of the exact solution
    function du =  Du(pt)
        x = pt(:,1); y = pt(:,2);
        du(:,1,1) = x*0;
        du(:,1,2) = x*0;
        du(:,2,1) = x*0;
        du(:,2,2) = x*0;
    end

    function dp =  Dp(pt)
        x = pt(:,1); y = pt(:,2);
        dp(:,1) = x*0;
        dp(:,2) = x*0;
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
        
        % note: zero
        uN(:,1) = x*0;
        uN(:,2) = x*0;
    end



    % g_RN = dudn + beta u
    function uR = g_R(pt)
        x = pt(:,1); y = pt(:,2);
        uR = 0*x + m0;
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

        % note: zero
        uN(:,1) = x*0;
        uN(:,2) = x*0;
    end

    function uN = g_Dn(pt)
        x = pt(:,1); y = pt(:,2);
        y1 = y + slope * x;

        u = exactu(pt);
        n = [-slope, -1];
        n = n / norm(n);

        % note: zero
        uN = x*0;
    end
    
end
