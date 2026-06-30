%% NONLINEARFOP2FEMRATEN3 convergence test for nonlinear FO P2 solver.
%
% Manufactured solution on the unit cube:
%
%   u = sin(pi*x)*sin(pi*y)*sin(pi*z),
%   v = cos(pi*x)*sin(pi*y)*sin(pi*z).
%
% This test uses Glen exponent n = 3, A = 1, and the regularized FO
% viscosity used by NonlinearFOP2.  The forcing is f = -div(stress(U)) for
% the exact nonlinear stress, evaluated analytically by the chain rule.

close all;

[node,elem] = cubemesh([0,1,0,1,0,1],0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');

epsReg = 1e-4;

pde.A = 1;
pde.n = 3;
pde.f = @(p) foforce(p,epsReg);
pde.g_D = @foexact;

option.printlevel = 0;
option.maxIt = 100;
option.tol = 1e-10;
option.residual_tol = 1e-10;
option.residual_check_threshold = 1e-4;
option.damping = 0.8;
option.eps_reg = epsReg;
option.quadorder = 5;

maxIt = 3;
h = zeros(maxIt,1);
N = zeros(maxIt,1);
errL2 = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
itStep = zeros(maxIt,1);
nlResidual = zeros(maxIt,1);

fprintf('\nConvergence table for NonlinearFOP2, n = 3\n');
fprintf('%10s  %10s  %12s  %7s  %12s  %7s  %4s  %10s\n',...
    '#Dof','h','L2 error','rate','H1 error','rate','it','nlres');
for k = 1:maxIt
    [soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);
    [errL2(k),errH1(k)] = foerror(node,elem,eqn.elem2dof,soln);
    h(k) = max(edgelength(node,elem));
    N(k) = length(soln.U);
    itStep(k) = info.itStep;
    nlResidual(k) = info.equationResidual(end);
    if k == 1
        rateL2 = NaN;
        rateH1 = NaN;
    else
        rateL2 = log(errL2(k-1)/errL2(k))/log(h(k-1)/h(k));
        rateH1 = log(errH1(k-1)/errH1(k))/log(h(k-1)/h(k));
    end
    fprintf('%10d  %10.3e  %12.4e  %7s  %12.4e  %7s  %4d  %10.3e\n',...
        N(k),h(k),errL2(k),ratestr(rateL2),errH1(k),...
        ratestr(rateH1),itStep(k),nlResidual(k));
    drawnow;

    if k < maxIt
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    end
end

function value = foexact(p)
x = p(:,1); y = p(:,2); z = p(:,3);
value = [sin(pi*x).*sin(pi*y).*sin(pi*z), ...
         cos(pi*x).*sin(pi*y).*sin(pi*z)];
end

function value = ratestr(rate)
if isnan(rate)
    value = '   --  ';
else
    value = sprintf('%7.2f',rate);
end
end

function value = fogradexact(p)
[d,~] = foderivatives(p);
value = [d.ux,d.uy,d.uz,d.vx,d.vy,d.vz];
end

function value = foforce(p,epsReg)
[d,dd] = foderivatives(p);

epsII = d.ux.^2+d.vy.^2+d.ux.*d.vy+...
    0.25*(d.uy+d.vx).^2+0.25*d.uz.^2+0.25*d.vz.^2;
strain = epsII+epsReg^2;
eta = 0.5*strain.^(-1/3);

Ex = epsderivative(d,dd,'x');
Ey = epsderivative(d,dd,'y');
Ez = epsderivative(d,dd,'z');
etaX = eta.*(-1/3).*Ex./strain;
etaY = eta.*(-1/3).*Ey./strain;
etaZ = eta.*(-1/3).*Ez./strain;

S1 = 4*d.ux+2*d.vy;
S2 = d.uy+d.vx;
S3 = d.uz;
T2 = 4*d.vy+2*d.ux;
T3 = d.vz;

S1x = 4*dd.uxx+2*dd.vxy;
S2x = dd.uxy+dd.vxx;
S2y = dd.uyy+dd.vxy;
S3z = dd.uzz;
T2y = 4*dd.vyy+2*dd.uxy;
T3z = dd.vzz;

div1 = etaX.*S1+eta.*S1x+etaY.*S2+eta.*S2y+...
    etaZ.*S3+eta.*S3z;
div2 = etaX.*S2+eta.*S2x+etaY.*T2+eta.*T2y+...
    etaZ.*T3+eta.*T3z;
value = [-div1,-div2];
end

function E = epsderivative(d,dd,coord)
switch coord
    case 'x'
        uxI = dd.uxx; uyI = dd.uxy; uzI = dd.uxz;
        vxI = dd.vxx; vyI = dd.vxy; vzI = dd.vxz;
    case 'y'
        uxI = dd.uxy; uyI = dd.uyy; uzI = dd.uyz;
        vxI = dd.vxy; vyI = dd.vyy; vzI = dd.vyz;
    case 'z'
        uxI = dd.uxz; uyI = dd.uyz; uzI = dd.uzz;
        vxI = dd.vxz; vyI = dd.vyz; vzI = dd.vzz;
end
E = 2*d.ux.*uxI+2*d.vy.*vyI+uxI.*d.vy+d.ux.*vyI+...
    0.5*(d.uy+d.vx).*(uyI+vxI)+0.5*d.uz.*uzI+...
    0.5*d.vz.*vzI;
end

function [d,dd] = foderivatives(p)
x = p(:,1); y = p(:,2); z = p(:,3);
sx = sin(pi*x); cx = cos(pi*x);
sy = sin(pi*y); cy = cos(pi*y);
sz = sin(pi*z); cz = cos(pi*z);
pi2 = pi^2;

d.ux = pi*cx.*sy.*sz;
d.uy = pi*sx.*cy.*sz;
d.uz = pi*sx.*sy.*cz;
d.vx = -pi*sx.*sy.*sz;
d.vy = pi*cx.*cy.*sz;
d.vz = pi*cx.*sy.*cz;

dd.uxx = -pi2*sx.*sy.*sz;
dd.uxy = pi2*cx.*cy.*sz;
dd.uxz = pi2*cx.*sy.*cz;
dd.uyy = -pi2*sx.*sy.*sz;
dd.uyz = pi2*sx.*cy.*cz;
dd.uzz = -pi2*sx.*sy.*sz;

dd.vxx = -pi2*cx.*sy.*sz;
dd.vxy = -pi2*sx.*cy.*sz;
dd.vxz = -pi2*sx.*sy.*cz;
dd.vyy = -pi2*cx.*sy.*sz;
dd.vyz = pi2*cx.*cy.*cz;
dd.vzz = -pi2*cx.*sy.*sz;
end

function [errL2,errH1] = foerror(node,elem,elem2dof,soln)
[Dlambda,volume] = gradbasis3(node,elem);
[lambda,weight] = quadpts3(5);
uh = soln.u;
vh = soln.v;
errL2sq = 0;
errH1sq = 0;

for q = 1:size(lambda,1)
    phi = p2basis3(lambda(q,:));
    Dphi = p2gradient3local(lambda(q,:),Dlambda);
    xq = lambda(q,1)*node(elem(:,1),:) + ...
         lambda(q,2)*node(elem(:,2),:) + ...
         lambda(q,3)*node(elem(:,3),:) + ...
         lambda(q,4)*node(elem(:,4),:);
    exact = foexact(xq);
    gradExact = fogradexact(xq);
    uhq = uh(elem2dof)*phi';
    vhq = vh(elem2dof)*phi';
    gradU = evalp2grad(uh,elem2dof,Dphi);
    gradV = evalp2grad(vh,elem2dof,Dphi);

    errL2sq = errL2sq + sum(weight(q)*volume.*...
        ((uhq-exact(:,1)).^2+(vhq-exact(:,2)).^2));
    du = gradU-gradExact(:,1:3);
    dv = gradV-gradExact(:,4:6);
    errH1sq = errH1sq + sum(weight(q)*volume.*...
        (sum(du.^2,2)+sum(dv.^2,2)));
end

errL2 = sqrt(errL2sq);
errH1 = sqrt(errH1sq);
end

function h = edgelength(node,elem)
edges = [elem(:,[1 2]); elem(:,[1 3]); elem(:,[1 4]); ...
         elem(:,[2 3]); elem(:,[2 4]); elem(:,[3 4])];
h = sqrt(sum((node(edges(:,1),:)-node(edges(:,2),:)).^2,2));
end

function phi = p2basis3(l)
phi = [l(1)*(2*l(1)-1), l(2)*(2*l(2)-1), ...
       l(3)*(2*l(3)-1), l(4)*(2*l(4)-1), ...
       4*l(1)*l(2), 4*l(1)*l(3), 4*l(1)*l(4), ...
       4*l(2)*l(3), 4*l(2)*l(4), 4*l(3)*l(4)];
end

function Dphi = p2gradient3local(l,Dlambda)
NT = size(Dlambda,1);
Dphi = zeros(NT,3,10);
Dphi(:,:,1) = (4*l(1)-1).*Dlambda(:,:,1);
Dphi(:,:,2) = (4*l(2)-1).*Dlambda(:,:,2);
Dphi(:,:,3) = (4*l(3)-1).*Dlambda(:,:,3);
Dphi(:,:,4) = (4*l(4)-1).*Dlambda(:,:,4);
Dphi(:,:,5) = 4*(l(1)*Dlambda(:,:,2)+l(2)*Dlambda(:,:,1));
Dphi(:,:,6) = 4*(l(1)*Dlambda(:,:,3)+l(3)*Dlambda(:,:,1));
Dphi(:,:,7) = 4*(l(1)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,1));
Dphi(:,:,8) = 4*(l(2)*Dlambda(:,:,3)+l(3)*Dlambda(:,:,2));
Dphi(:,:,9) = 4*(l(2)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,2));
Dphi(:,:,10) = 4*(l(3)*Dlambda(:,:,4)+l(4)*Dlambda(:,:,3));
end

function gradValue = evalp2grad(value,elem2dof,Dphi)
NT = size(Dphi,1);
localValue = reshape(value(elem2dof),NT,1,10);
gradValue = zeros(NT,3);
for a = 1:10
    gradValue = gradValue+Dphi(:,:,a).*localValue(:,:,a);
end
end
