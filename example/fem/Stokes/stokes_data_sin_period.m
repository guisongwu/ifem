function pde = stokes_data_sin_period
%% STOKESDATA2 data for Stokes equations
%
% The solution u is a polynomial satisfying the zero Dirichlet boundary
% condition. 
%
% Created by Ming Wang.
%
    pde = struct('f', @f, ...
                 'fp', @fp, ...
                 'exactp', @exactp, ...
                 'exactu',@exactu, ...
                 'exactux',@exactux, ...
                 'exactuy',@exactuy, ...
                 'g_D', @exactu, 'g_N', @g_N, 'g_R',@g_R, 'g_RN', @g_RN, 'g_Dn', @g_Dn);
%

global slope dbg_case h;

Power = @(x,n) power(x,n);
Sin = @(x) sin(x);
Cos = @(x) cos(x);
Pi = pi;

if dbg_case == 2
    xshift = 3*h;
else
    xshift = 0;
end

% (1 + slope*x + y)*Sin(2*Pi*x)
% (1 - Power(slope*x + y,2)/2.)*Cos(2*Pi*x)
% (1 - slope*x - y)*Cos(2*Pi*x)
% -(slope*Cos(2*Pi*x)) - 4*Pi*slope*Cos(2*Pi*x) - 2*Pi*(1 - slope*x - y).*Sin(2*Pi*x) + 4*Power(Pi,2)*(1 + slope*x + y).*Sin(2*Pi*x)
% Power(slope,2)*Cos(2*Pi*x) + 4*Power(Pi,2)*(1 - Power(slope*x + y,2)/2.)*Cos(2*Pi*x) - 4*Pi*slope*(slope*x + y).*Sin(2*Pi*x)
% (-(slope*x) - y)*Cos(2*Pi*x) + 2*Pi*(1 + slope*x + y)*Cos(2*Pi*x) + slope.*Sin(2*Pi*x)

%% subfunction
    function z = exactux(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z = (1 + slope*x + y).*Sin(2*Pi*x);
    end
    function z = exactuy(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z = (1 - Power(slope*x + y,2)/2.).*Cos(2*Pi*x);
    end

    function z = exactu(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z(:,1) = exactux(pt);
        z(:,2) = exactuy(pt);
    end

    function z = exactp(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z = (1 - slope*x - y).*Cos(2*Pi*x);
    end

    function z = f(pt) % load data (right hand side function)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z(:,1) = -(slope.*Cos(2*Pi*x)) - 4*Pi*slope.*Cos(2*Pi*x) - 2*Pi*(1 - slope*x - y).*Sin(2*Pi*x) + 4*Power(Pi,2)*(1 + slope*x + y).*Sin(2*Pi*x);
        z(:,2) = Power(slope,2).*Cos(2*Pi*x) + 4*Power(Pi,2)*(1 - Power(slope*x + y,2)/2.).*Cos(2*Pi*x) - 4*Pi*slope*(slope*x + y).*Sin(2*Pi*x);
    end

    function z = fp(pt) % load data (right hand side function)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        z = (-(slope*x) - y).*Cos(2*Pi*x) + 2*Pi*(1 + slope*x + y).*Cos(2*Pi*x) + slope.*Sin(2*Pi*x);
    end



    
    
    
    
% 2*Pi*(1 + slope*x + y).*Cos(2*Pi*x) + slope.*Sin(2*Pi*x)
% Sin(2*Pi*x)
% -(slope*(slope*x + y).*Cos(2*Pi*x)) - 2*Pi*(1 - Power(slope*x + y,2)/2.).*Sin(2*Pi*x)
% (-(slope*x) - y).*Cos(2*Pi*x)
% -(slope.*Cos(2*Pi*x)) - 2*Pi*(1 - slope*x - y).*Sin(2*Pi*x)
% -Cos(2*Pi*x)
   
    % Derivative of the exact solution
    function du =  Du(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        du(:,1,1) = 2*Pi*(1 + slope*x + y).*Cos(2*Pi*x) + slope.*Sin(2*Pi*x);
        du(:,1,2) = Sin(2*Pi*x);
        du(:,2,1) = -(slope*(slope*x + y).*Cos(2*Pi*x)) - 2*Pi*(1 - Power(slope*x + y,2)/2.).*Sin(2*Pi*x);
        du(:,2,2) = (-(slope*x) - y).*Cos(2*Pi*x);
    end

    function dp =  Dp(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

        dp(:,1) = -(slope.*Cos(2*Pi*x)) - 2*Pi*(1 - slope*x - y).*Sin(2*Pi*x);
        dp(:,2) = -Cos(2*Pi*x);
    end




    % g_N = dudn
    function uN = g_N(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

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
        x = x - xshift;
        y = y - slope * xshift;

        uR = 0*x + 1 + 0.5*Cos(2*Pi*x);
    end
    
    function uN = g_RN(pt)
        x = pt(:,1); y = pt(:,2);
        x = x - xshift;
        y = y - slope * xshift;

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
        x = x - xshift;
        y = y - slope * xshift;

        y1 = y + slope * x;

        u = exactu(pt);
        n = [-slope, -1];
        n = n / norm(n);

        uN = u(:,1) * n(1) + u(:,2) * n(2);
    end
    
end
