function pde = NonlinearStokesMMSData(L,H,slope,eps_reg)
%% NONLINEARSTOKESMMSDATA Manufactured solution for nonlinear Stokes.
%
% The physical coordinates are (x,z), and q = z+slope*x.  The stream
% function psi = sin(2*pi*x/L)*q^3*(H-q)^3 gives a divergence-free,
% periodic velocity that vanishes at the bed q=0.

if nargin < 4, eps_reg = 1e-4; end

pde.A = 1;
pde.n = 3;
pde.beta = 10;
pde.m = 1/3;
pde.f = @bodyforce;
pde.g_N = @toptraction;
pde.exactu = @exactu;
pde.exactux = @exactux;
pde.exactuz = @exactuz;
pde.exactp = @exactp;

k = 2*pi/L;
nGlen = pde.n;
A = pde.A;

    function u = exactu(pt)
        [ux,uz] = velocity(pt);
        u = [ux,uz];
    end

    function ux = exactux(pt)
        [ux,~] = velocity(pt);
    end

    function uz = exactuz(pt)
        [~,uz] = velocity(pt);
    end

    function p = exactp(pt)
        x = pt(:,1);
        q = pt(:,2)+slope*x;
        p = cos(k*x).*(1+q);
    end

    function f = bodyforce(pt)
        d = fields(pt);
        exponent = (1-nGlen)/(2*nGlen);
        eta = 0.5*A^(-1/nGlen)*d.Q.^exponent;
        etax = eta.*exponent.*d.Qx./d.Q;
        etaz = eta.*exponent.*d.Qz./d.Q;

        fx = d.px-2*(etax.*d.exx+eta.*d.exxx+...
                     etaz.*d.exz+eta.*d.exzz);
        fz = d.pz-2*(etax.*d.exz+eta.*d.exzx+...
                     etaz.*d.ezz+eta.*d.ezzz);
        f = [fx,fz];
    end

    function traction = toptraction(pt)
        d = fields(pt);
        exponent = (1-nGlen)/(2*nGlen);
        eta = 0.5*A^(-1/nGlen)*d.Q.^exponent;
        sigma11 = 2*eta.*d.exx-d.p;
        sigma12 = 2*eta.*d.exz;
        sigma22 = 2*eta.*d.ezz-d.p;
        normal = [slope,1]/sqrt(1+slope^2);
        traction = [sigma11*normal(1)+sigma12*normal(2),...
                    sigma12*normal(1)+sigma22*normal(2)];
    end

    function [ux,uz] = velocity(pt)
        x = pt(:,1);
        q = pt(:,2)+slope*x;
        [F,F1] = bedpolynomial(q);
        S = sin(k*x);
        C = cos(k*x);
        ux = S.*F1;
        uz = -k*C.*F-slope*S.*F1;
    end

    function d = fields(pt)
        x = pt(:,1);
        q = pt(:,2)+slope*x;
        [F,F1,F2,F3] = bedpolynomial(q);
        S = sin(k*x);
        C = cos(k*x);
        s = slope;

        ax = k*C.*F1+s*S.*F2;
        az = S.*F2;
        bx = k^2*S.*F-2*s*k*C.*F1-s^2*S.*F2;
        bz = -k*C.*F1-s*S.*F2;

        axx = -k^2*S.*F1+2*s*k*C.*F2+s^2*S.*F3;
        axz = k*C.*F2+s*S.*F3;
        azz = S.*F3;
        bxx = k^3*C.*F+3*s*k^2*S.*F1-...
              3*s^2*k*C.*F2-s^3*S.*F3;
        bxz = k^2*S.*F1-2*s*k*C.*F2-s^2*S.*F3;
        bzz = -k*C.*F2-s*S.*F3;

        d.exx = ax;
        d.ezz = bz;
        d.exz = 0.5*(az+bx);
        d.exxx = axx;
        d.exxz = axz;
        d.ezzx = bxz;
        d.ezzz = bzz;
        d.exzx = 0.5*(axz+bxx);
        d.exzz = 0.5*(azz+bxz);

        epsII = d.exx.^2+d.exz.^2;
        d.Q = epsII+eps_reg^2;
        d.Qx = 2*d.exx.*d.exxx+2*d.exz.*d.exzx;
        d.Qz = 2*d.exx.*d.exxz+2*d.exz.*d.exzz;

        d.p = C.*(1+q);
        d.px = -k*S.*(1+q)+s*C;
        d.pz = C;
    end

    function [F,F1,F2,F3] = bedpolynomial(q)
        F = H^3*q.^3-3*H^2*q.^4+3*H*q.^5-q.^6;
        F1 = 3*H^3*q.^2-12*H^2*q.^3+15*H*q.^4-6*q.^5;
        if nargout > 2
            F2 = 6*H^3*q-36*H^2*q.^2+60*H*q.^3-30*q.^4;
            F3 = 6*H^3-72*H^2*q+180*H*q.^2-120*q.^3;
        end
    end
end
