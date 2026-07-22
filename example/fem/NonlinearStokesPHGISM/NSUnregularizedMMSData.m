function pde = NSUnregularizedMMSData(L,H,slope)
%% NONLINEARSTOKESUNREGULARIZEDMMSDATA Exact shear solution without eps_reg.
%
% The exact solution is a tangential shear flow in the slab coordinates
% q = z+slope*x:
%
%   u_x = U0+gamma*q,   u_z = -slope*(U0+gamma*q),   p = 0.
%
% The strain rate is constant and nonzero, so the unregularized Glen
% viscosity is finite.  The bed sliding coefficient beta is chosen so the
% Robin sliding boundary condition is satisfied exactly.

if nargin < 1, L = 1; end
if nargin < 2, H = 1; end
if nargin < 3, slope = 0.1; end

pde.A = 1;
pde.n = 3;
pde.m = 1/3;
pde.f = @bodyforce;
pde.g_N = @toptraction;
pde.beta = @bedbeta;
pde.exactu = @exactu;
pde.exactux = @exactux;
pde.exactuz = @exactuz;
pde.exactp = @exactp;

U0 = 1+0*L;
gamma = 1+0*H;

tangent = [1,-slope]/sqrt(1+slope^2);
topNormal = [slope,1]/sqrt(1+slope^2);
bedNormal = [-slope,-1]/sqrt(1+slope^2);
sigma = stressmatrix;
bedTangentialTraction = tangent*(sigma*bedNormal');
bedTangentialVelocity = U0*sqrt(1+slope^2);
betaValue = -bedTangentialTraction ...
    /(abs(bedTangentialVelocity)^(pde.m-1)*bedTangentialVelocity);

if betaValue <= 0
    error('iFEM:UnregularizedMMSBeta',...
        'The manufactured bed beta must be positive.');
end

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
        p = zeros(size(pt,1),1);
    end

    function f = bodyforce(pt)
        f = zeros(size(pt,1),2);
    end

    function traction = toptraction(pt)
        nt = size(pt,1);
        tractionVector = sigma*topNormal';
        traction = repmat(tractionVector',nt,1);
    end

    function beta = bedbeta(pt)
        beta = betaValue+zeros(size(pt,1),1);
    end

    function [ux,uz] = velocity(pt)
        q = pt(:,2)+slope*pt(:,1);
        speed = U0+gamma*q;
        ux = speed;
        uz = -slope*speed;
    end

    function sigmaMatrix = stressmatrix
        ux_x = gamma*slope;
        ux_z = gamma;
        uz_x = -slope^2*gamma;
        uz_z = -slope*gamma;

        exx = ux_x;
        ezz = uz_z;
        exz = 0.5*(ux_z+uz_x);
        epsII = 0.5*(exx^2+ezz^2+2*exz^2);
        eta = 0.5*pde.A^(-1/pde.n)*epsII^((1-pde.n)/(2*pde.n));

        sigmaMatrix = [2*eta*exx, 2*eta*exz;...
                       2*eta*exz, 2*eta*ezz];
    end
end
