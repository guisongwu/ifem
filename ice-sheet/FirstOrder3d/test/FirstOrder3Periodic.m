function [soln,eqn,info,node,elem,bdFlag] = FirstOrder3Periodic(option)
%% FIRSTORDER3PERIODIC FO P2 solver on a periodic rectangular box.
%
% Solve the nonlinear three-dimensional FO/BP model on
%
%   [0,4] x [0,4] x [0,1].
%
% Boundary conditions:
%   top z = 1: zero horizontal traction;
%   front/back y = 0,4: periodic;
%   left/right x = 0,4: periodic;
%   bottom z = 0: zero normal penetration is built into the FO horizontal
%                 velocity formulation, with tangential linear sliding m=1.
%
% Example:
%
%   [soln,eqn,info,node,elem,bdFlag] = FirstOrder3Periodic;

if nargin < 1, option = struct; end

box = [0,4,0,4,0,1];
h = getoption(option,'h',[1,1,0.5]);
if isscalar(h)
    h = [h,h,min(h,box(6)-box(5))];
end
[node,elem] = cubemesh(box,h);

bdFlag = setboundary3(node,elem,...
    'Neumann','z==1',...
    'Robin','z==0');
if ~any(bdFlag(:)==2) || ~any(bdFlag(:)==3)
    error('iFEM:FOBoundary',...
        'The mesh must contain both z=1 top faces and z=0 bottom faces.');
end

pde.A = getoption(option,'A',1);
pde.n = getoption(option,'n',3);
pde.rho = getoption(option,'rho',910);
pde.gravity = getoption(option,'gravity',9.81);
pde.beta = getoption(option,'beta',1e3);
pde.m = 1;
pde.gradS = getoption(option,'gradS',@(p) constantgradient(p,[1e-2,0]));
pde.g_N = @(p) zeros(size(p,1),2);

option = setoption(option,'bed_condition','sliding');
option = setoption(option,'periodic',[1 2]);
option = setoption(option,'periodicBox',[0 4; 0 4; 0 1]);
option = setoption(option,'printlevel',1);
% option = setoption(option,'maxIt',50);
option = setoption(option,'maxIt',100);
option = setoption(option,'tol',1e-8);
option = setoption(option,'residual_tol',1e-8);
option = setoption(option,'damping',0.7);

[soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);

if nargout == 0
    speed = sqrt(soln.u(1:size(node,1)).^2+soln.v(1:size(node,1)).^2);
    fprintf('FO rectangular periodic solve: %d dofs, %d iterations, max speed %.4e\n',...
        length(soln.U),info.itStep,max(speed));
    clear soln eqn info node elem bdFlag
end

end

function value = constantgradient(p,g)
value = repmat(g,size(p,1),1);
end

function value = getoption(option,name,defaultValue)
if isfield(option,name)
    value = option.(name);
else
    value = defaultValue;
end
end

function option = setoption(option,name,value)
if ~isfield(option,name), option.(name) = value; end
end
