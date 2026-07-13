%% FIRSTORDER3VISUALIZE Visualize a manufactured 3-D FO P2 solution.
%
% The script solves the same manufactured problem used by
% FirstOrder3ConverRate on a moderately refined cube, then plots boundary
% speed, a middle z-slice, and horizontal velocity vectors.

close all;

exportPng = false;
moduleDir = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(moduleDir,'output',mfilename);

[node,elem] = cubemesh([0,1,0,1,0,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet');
for k = 1:2
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end

pde.A = 1;
pde.n = 1;
pde.f = @foforce;
pde.g_D = @foexact;

option.printlevel = 1;
option.maxIt = 30;
option.tol = 1e-12;
option.residual_tol = 1e-12;
option.damping = 1;
option.quadorder = 5;

[soln,eqn,info] = NonlinearFOP2(node,elem,bdFlag,pde,option);
fprintf('Visualization solve: dof=%d, it=%d, residual=%.3e\n',...
    length(soln.U),info.itStep,info.equationResidual(end));

N = size(node,1);
udofNode = [node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2];
speedDof = sqrt(soln.u.^2+soln.v.^2);
speedNode = speedDof(1:N);

fig1 = figure('Name','boundary speed');
[~,bdFace] = findboundary3(elem);
trisurf(bdFace,node(:,1),node(:,2),node(:,3),speedNode,...
    'EdgeColor',[0.65 0.65 0.65],'FaceColor','interp');
axis equal tight;
view(35,24);
colorbar;
xlabel('x'); ylabel('y'); zlabel('z');
title('|U_h| on boundary');
exportfigure(fig1,outputDir,'FirstOrder3-boundary-speed.png',exportPng);

ns = 41;
[xs,ys] = meshgrid(linspace(0,1,ns),linspace(0,1,ns));
zs = 0.5*ones(size(xs));
Fu = scatteredInterpolant(udofNode(:,1),udofNode(:,2),udofNode(:,3),...
    soln.u,'natural','none');
Fv = scatteredInterpolant(udofNode(:,1),udofNode(:,2),udofNode(:,3),...
    soln.v,'natural','none');
us = Fu(xs,ys,zs);
vs = Fv(xs,ys,zs);
speeds = sqrt(us.^2+vs.^2);

fig2 = figure('Name','mid-slice speed');
surf(xs,ys,zs,speeds,'EdgeColor','none');
axis equal tight;
view(2);
colorbar;
xlabel('x'); ylabel('y');
title('|U_h| on z = 0.5');

hold on;
stride = 4;
quiver3(xs(1:stride:end,1:stride:end),...
        ys(1:stride:end,1:stride:end),...
        zs(1:stride:end,1:stride:end),...
        us(1:stride:end,1:stride:end),...
        vs(1:stride:end,1:stride:end),...
        zeros(size(us(1:stride:end,1:stride:end))),...
        'k');
hold off;
exportfigure(fig2,outputDir,'FirstOrder3-mid-slice-speed.png',exportPng);

fig3 = figure('Name','exact-vs-numerical slice error');
exactSlice = foexact([xs(:),ys(:),zs(:)]);
exactSpeed = reshape(sqrt(exactSlice(:,1).^2+exactSlice(:,2).^2),ns,ns);
surf(xs,ys,zs,abs(speeds-exactSpeed),'EdgeColor','none');
axis equal tight;
view(2);
colorbar;
xlabel('x'); ylabel('y');
title('speed error on z = 0.5');
exportfigure(fig3,outputDir,'FirstOrder3-mid-slice-error.png',exportPng);

function value = foexact(p)
x = p(:,1); y = p(:,2); z = p(:,3);
value = [sin(pi*x).*sin(pi*y).*sin(pi*z), ...
         cos(pi*x).*sin(pi*y).*sin(pi*z)];
end

function value = foforce(p)
x = p(:,1); y = p(:,2); z = p(:,3);
f1 = pi^2*sin(pi*x).*sin(pi*z).*...
    (3*sin(pi*y)+1.5*cos(pi*y));
f2 = pi^2*cos(pi*x).*sin(pi*z).*...
    (3*sin(pi*y)-1.5*cos(pi*y));
value = [f1,f2];
end

function exportfigure(fig,outputDir,fileName,exportPng)
if ~exportPng, return; end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
exportgraphics(fig,fullfile(outputDir,fileName),'Resolution',200);
end
