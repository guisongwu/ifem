%% NONLINEARFOP2FEMRATE convergence test for the FO P2 solver.
%
% Manufactured solution on the unit cube:
%
%   u = sin(pi*x)*sin(pi*y)*sin(pi*z),
%   v = cos(pi*x)*sin(pi*y)*sin(pi*z).
%
% With n = 1 and A = 1, Glen viscosity is eta = 1/2 and the FO operator is
% linear.  P2 elements should show about third order in L2 and second order
% in the H1 seminorm on uniform refinements.

close all;

[node,elem] = cubemesh([0,1,0,1,0,1],0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');

pde.A = 1;
pde.n = 1;
pde.f = @foforce;
pde.g_D = @foexact;

option.printlevel = 0;
option.maxIt = 30;
option.tol = 1e-12;
option.residual_tol = 1e-12;
option.damping = 1;
option.quadorder = 5;

maxIt = 4;
h = zeros(maxIt,1);
N = zeros(maxIt,1);
errL2 = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
itStep = zeros(maxIt,1);

fprintf('\nConvergence table for NonlinearFOP2\n');
fprintf('%10s  %10s  %12s  %7s  %12s  %7s  %4s\n',...
    '#Dof','h','L2 error','rate','H1 error','rate','it');
for k = 1:maxIt
    [soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);
    [errL2(k),errH1(k)] = foerror(node,elem,eqn.elem2dof,soln);
    h(k) = max(edgelength(node,elem));
    N(k) = length(soln.U);
    itStep(k) = info.itStep;
    if k == 1
        rateL2 = NaN;
        rateH1 = NaN;
    else
        rateL2 = log(errL2(k-1)/errL2(k))/log(h(k-1)/h(k));
        rateH1 = log(errH1(k-1)/errH1(k))/log(h(k-1)/h(k));
    end
    fprintf('%10d  %10.3e  %12.4e  %7s  %12.4e  %7s  %4d\n',...
        N(k),h(k),errL2(k),ratestr(rateL2),errH1(k),...
        ratestr(rateH1),itStep(k));
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
x = p(:,1); y = p(:,2); z = p(:,3);
du = [pi*cos(pi*x).*sin(pi*y).*sin(pi*z), ...
      pi*sin(pi*x).*cos(pi*y).*sin(pi*z), ...
      pi*sin(pi*x).*sin(pi*y).*cos(pi*z)];
dv = [-pi*sin(pi*x).*sin(pi*y).*sin(pi*z), ...
       pi*cos(pi*x).*cos(pi*y).*sin(pi*z), ...
       pi*cos(pi*x).*sin(pi*y).*cos(pi*z)];
value = [du,dv];
end

function value = foforce(p)
x = p(:,1); y = p(:,2); z = p(:,3);
f1 = pi^2*sin(pi*x).*sin(pi*z).*...
    (3*sin(pi*y)+1.5*cos(pi*y));
f2 = pi^2*cos(pi*x).*sin(pi*z).*...
    (3*sin(pi*y)-1.5*cos(pi*y));
value = [f1,f2];
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
