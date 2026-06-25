function pde = NonlinearStokesBMMSData(L,H,bedAmplitude,slope,eps_reg)
%% NONLINEARSTOKESBMMSDATA Manufactured solution on an ISMIP-HOM B bed.
%
% The geometry is the flowline ISMIP-HOM B geometry
%
%   surface(x) = -slope*x,
%   bed(x)     = surface(x)-H+bedAmplitude*sin(2*pi*x/L).
%
% The stream function is periodic in x and uses the mapped thickness
% coordinate r=(z-bed)/(surface-bed).  Since psi is constant on the bed,
% the exact velocity is impermeable there.  Its tangential component is
% nonzero, so the basal Robin coefficient is recovered from the exact
% stress and the sliding condition.

if nargin < 5, eps_reg = 1e-2; end

pde.A = 1;
pde.n = 3;
pde.m = 1;
pde.beta = @bedbeta;
pde.f = @bodyforce;
pde.g_N = @toptraction;
pde.exactu = @exactu;
pde.exactux = @exactux;
pde.exactuz = @exactuz;
pde.exactp = @exactp;
pde.geometry.surface = @surface;
pde.geometry.bed = @bed;

k = 2*pi/L;
nGlen = pde.n;
A = pde.A;
streamAmplitude = 1;
waveAmplitude = 0.25;
shapeAmplitude = 1;
diffStep = 1e-5*max(1,min(L,H));

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
        g = geometry(pt);
        p = cos(k*g.x).*(1+g.r);
    end

    function f = bodyforce(pt)
        h = diffStep;
        [s11px,s12px,~] = stress(pt+[h,0]);
        [s11mx,s12mx,~] = stress(pt-[h,0]);
        [~,s12pz,s22pz] = stress(pt+[0,h]);
        [~,s12mz,s22mz] = stress(pt-[0,h]);

        divx = (s11px-s11mx)/(2*h)+(s12pz-s12mz)/(2*h);
        divz = (s12px-s12mx)/(2*h)+(s22pz-s22mz)/(2*h);
        f = [-divx,-divz];
    end

    function traction = toptraction(pt)
        [s11,s12,s22] = stress(pt);
        n = topnormal(pt);
        traction = [s11.*n(:,1)+s12.*n(:,2),...
                    s12.*n(:,1)+s22.*n(:,2)];
    end

    function beta = bedbeta(pt)
        [s11,s12,s22] = stress(pt);
        t = bedtangent(pt);
        n = [t(:,2),-t(:,1)];
        [ux,uz] = velocity(pt);
        ut = ux.*t(:,1)+uz.*t(:,2);
        tau = t(:,1).*(s11.*n(:,1)+s12.*n(:,2)) + ...
              t(:,2).*(s12.*n(:,1)+s22.*n(:,2));
        beta = -tau./ut;
    end

    function [s11,s12,s22] = stress(pt)
        [uxx,uxz,uzx,uzz] = velocitygradient(pt);
        exx = uxx;
        ezz = uzz;
        exz = 0.5*(uxz+uzx);
        epsII = exx.^2+exz.^2;
        eta = 0.5*A^(-1/nGlen)*(epsII+eps_reg^2).^...
            ((1-nGlen)/(2*nGlen));
        p = exactp(pt);
        s11 = 2*eta.*exx-p;
        s12 = 2*eta.*exz;
        s22 = 2*eta.*ezz-p;
    end

    function [uxx,uxz,uzx,uzz] = velocitygradient(pt)
        h = diffStep;
        [uxpx,uzpx] = velocity(pt+[h,0]);
        [uxmx,uzmx] = velocity(pt-[h,0]);
        [uxpz,uzpz] = velocity(pt+[0,h]);
        [uxmz,uzmz] = velocity(pt-[0,h]);
        uxx = (uxpx-uxmx)/(2*h);
        uzx = (uzpx-uzmx)/(2*h);
        uxz = (uxpz-uxmz)/(2*h);
        uzz = (uzpz-uzmz)/(2*h);
    end

    function [ux,uz] = velocity(pt)
        g = geometry(pt);
        [G,G1] = shapefunction(g.r);
        P = 1+waveAmplitude*sin(k*g.x);
        P1 = waveAmplitude*k*cos(k*g.x);
        ux = streamAmplitude*P.*G1.*g.rz;
        uz = -streamAmplitude*(P1.*G+P.*G1.*g.rx);
    end

    function g = geometry(pt)
        x = pt(:,1);
        z = pt(:,2);
        s = surfacex(x);
        b = bedx(x);
        T = s-b;
        sp = -slope+0*x;
        bp = -slope+bedAmplitude*k*cos(k*x);
        Tp = sp-bp;
        r = (z-b)./T;
        g.x = x;
        g.z = z;
        g.r = r;
        g.T = T;
        g.rx = -(bp+r.*Tp)./T;
        g.rz = 1./T;
    end

    function s = surface(pt)
        s = surfacex(pt(:,1));
    end

    function b = bed(pt)
        b = bedx(pt(:,1));
    end

    function s = surfacex(x)
        s = -slope*x;
    end

    function b = bedx(x)
        b = surfacex(x)-H+bedAmplitude*sin(k*x);
    end

    function t = bedtangent(pt)
        x = pt(:,1);
        bp = -slope+bedAmplitude*k*cos(k*x);
        t = [ones(size(x)),bp];
        t = t./sqrt(sum(t.^2,2));
    end

    function n = topnormal(pt)
        x = pt(:,1);
        sp = -slope+0*x;
        n = [-sp,ones(size(x))];
        n = n./sqrt(sum(n.^2,2));
    end

    function [G,G1] = shapefunction(r)
        G = r+shapeAmplitude*r.^2.*(1-r).^2;
        G1 = 1+shapeAmplitude*(2*r-6*r.^2+4*r.^3);
    end
end
