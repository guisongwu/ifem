function pde = NSMMSData(L,W,H,slope,eps_reg)
%% NSMMSDATA Manufactured solution for 3-D nonlinear Stokes.
%
% The solution is periodic in x and y and divergence-free.  With
% q = z+slope*x and psi = sin(2*pi*x/L)*sin(2*pi*y/W)*q^2*sin(pi*q/H),
%
%   u = curl(0,psi,0) = (-psi_z,0,psi_x).
%
% This gives a genuinely three-dimensional velocity field.  The body force
% is generated from the same regularized Glen viscosity as the solver.  The
% velocity vanishes at q = 0, so the manufactured test uses a no-slip bed.

if nargin < 5, eps_reg = 1e-4; end
if nargin < 4, slope = 0; end
if nargin < 2 || isempty(W), W = L; end

pde.A = 1;
pde.n = 3;
pde.beta = 0;
pde.m = 1;
pde.eps_reg = eps_reg;
pde.f = @bodyforce;
pde.g_N = @toptraction;
pde.exactu = @exactu;
pde.exactux = @exactux;
pde.exactuy = @exactuy;
pde.exactuz = @exactuz;
pde.exactp = @exactp;

kx = 2*pi/L;
ky = 2*pi/W;
gammaP = 0.3;
fdStep = 1e-5*max([L,W,H,1]);

    function u = exactu(pt)
        [ux,uy,uz] = velocity(pt);
        u = [ux,uy,uz];
    end

    function ux = exactux(pt)
        [ux,~,~] = velocity(pt);
    end

    function uy = exactuy(pt)
        uy = zeros(size(pt,1),1);
    end

    function uz = exactuz(pt)
        [~,~,uz] = velocity(pt);
    end

    function p = exactp(pt)
        [~,Cx,~,Cy,q] = trigfields(pt);
        p = Cx.*Cy.*pressureprofile(q);
    end

    function f = bodyforce(pt)
        divTau = viscousstressdivergence(pt);
        gradP = pressuregradient(pt);
        f = -divTau+gradP;
    end

    function traction = toptraction(pt)
        tau = viscousstress(pt);
        p = exactp(pt);
        normal = [slope,0,1]/sqrt(1+slope^2);
        traction = [(tau(:,1)-p)*normal(1)+tau(:,2)*normal(2)+...
                    tau(:,3)*normal(3),...
                    tau(:,2)*normal(1)+(tau(:,4)-p)*normal(2)+...
                    tau(:,5)*normal(3),...
                    tau(:,3)*normal(1)+tau(:,5)*normal(2)+...
                    (tau(:,6)-p)*normal(3)];
    end

    function [ux,uy,uz] = velocity(pt)
        [Sx,Cx,Sy,~,q] = trigfields(pt);
        [F,F1] = bedpolynomial(q);
        ux = -Sx.*Sy.*F1;
        uy = zeros(size(ux));
        uz = kx*Cx.*Sy.*F+slope*Sx.*Sy.*F1;
    end

    function g = gradientfields(pt)
        [Sx,Cx,Sy,Cy,q] = trigfields(pt);
        [F,F1,F2] = bedpolynomial(q);
        s = slope;
        g.uxx = -kx*Cx.*Sy.*F1-s*Sx.*Sy.*F2;
        g.uxy = -ky*Sx.*Cy.*F1;
        g.uxz = -Sx.*Sy.*F2;
        g.uzx = -kx^2*Sx.*Sy.*F+2*s*kx*Cx.*Sy.*F1+...
            s^2*Sx.*Sy.*F2;
        g.uzy = ky*(kx*Cx.*Cy.*F+s*Sx.*Cy.*F1);
        g.uzz = kx*Cx.*Sy.*F1+s*Sx.*Sy.*F2;
    end

    function tau = viscousstress(pt)
        g = gradientfields(pt);
        exx = g.uxx;
        eyy = zeros(size(exx));
        ezz = g.uzz;
        exy = 0.5*g.uxy;
        exz = 0.5*(g.uxz+g.uzx);
        eyz = 0.5*g.uzy;
        epsII = 0.5*(exx.^2+eyy.^2+ezz.^2+...
            2*(exy.^2+exz.^2+eyz.^2));
        etaGlen = 0.5*pde.A^(-1/pde.n).*...
            (epsII+eps_reg^2).^((1-pde.n)/(2*pde.n));
        tau = [2*etaGlen.*exx,2*etaGlen.*exy,2*etaGlen.*exz,...
               2*etaGlen.*eyy,2*etaGlen.*eyz,2*etaGlen.*ezz];
    end

    function divTau = viscousstressdivergence(pt)
        dTauDx = fourthorderderivative(pt,1);
        dTauDy = fourthorderderivative(pt,2);
        dTauDz = fourthorderderivative(pt,3);
        divTau = [dTauDx(:,1)+dTauDy(:,2)+dTauDz(:,3),...
                  dTauDx(:,2)+dTauDy(:,4)+dTauDz(:,5),...
                  dTauDx(:,3)+dTauDy(:,5)+dTauDz(:,6)];
    end

    function dTau = fourthorderderivative(pt,dim)
        e = zeros(1,3);
        e(dim) = fdStep;
        dTau = (-viscousstress(pt+2*e)+8*viscousstress(pt+e)-...
                8*viscousstress(pt-e)+viscousstress(pt-2*e))/(12*fdStep);
    end

    function gradP = pressuregradient(pt)
        [Sx,Cx,Sy,Cy,q] = trigfields(pt);
        P = pressureprofile(q);
        P1 = pressureprofilederivative(q);
        gradP = [-kx*Sx.*Cy.*P+slope*Cx.*Cy.*P1,...
                 -ky*Cx.*Sy.*P,...
                 Cx.*Cy.*P1];
    end

    function [Sx,Cx,Sy,Cy,q] = trigfields(pt)
        x = pt(:,1);
        y = pt(:,2);
        q = pt(:,3)+slope*x;
        Sx = sin(kx*x);
        Cx = cos(kx*x);
        Sy = sin(ky*y);
        Cy = cos(ky*y);
    end

    function [F,F1,F2,F3] = bedpolynomial(q)
        a = pi/H;
        F = q.^2.*sin(a*q);
        F1 = 2*q.*sin(a*q)+a*q.^2.*cos(a*q);
        if nargout > 2
            F2 = 2*sin(a*q)+4*a*q.*cos(a*q)-a^2*q.^2.*sin(a*q);
            F3 = 6*a*cos(a*q)-6*a^2*q.*sin(a*q)-...
                a^3*q.^2.*cos(a*q);
        end
    end

    function P = pressureprofile(q)
        P = 1+q+gammaP*sin(pi*q/H);
    end

    function P1 = pressureprofilederivative(q)
        P1 = 1+gammaP*(pi/H)*cos(pi*q/H);
    end
end
