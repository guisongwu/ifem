function figs = FirstOrder3PeriodicVisualize(option,soln,eqn,info,node,elem,bdFlag)
%% FIRSTORDER3PERIODICVISUALIZE Visualize the periodic-box FO solution.
%
% Example:
%
%   FirstOrder3PeriodicVisualize
%
% Reuse an existing solve:
%
%   [soln,eqn,info,node,elem,bdFlag] = FirstOrder3Periodic(option);
%   FirstOrder3PeriodicVisualize(option,soln,eqn,info,node,elem,bdFlag)
%
% To save figures:
%
%   option.exportPng = true;
%   FirstOrder3PeriodicVisualize(option)

if nargin < 1, option = struct; end

exportPng = getoption(option,'exportPng',false);
moduleDir = fileparts(fileparts(mfilename('fullpath')));
outputDir = getoption(option,'outputDir',...
    fullfile(moduleDir,'output',mfilename));
visible = getoption(option,'visible','on');

if nargin < 2 || isempty(soln)
    solveOption = option;
    solveOption.printlevel = getoption(option,'printlevel',1);
    solveOption.h = getoption(option,'h',2);
    [soln,eqn,info,node,elem,bdFlag] = FirstOrder3Periodic(solveOption);
elseif nargin < 7
    error('iFEM:FOVisualizationInput',...
        'Provide soln, eqn, info, node, elem, and bdFlag together.');
end

fprintf('Rectangular FO visualization: dof=%d, it=%d, residual=%.3e\n',...
    length(soln.U),info.itStep,info.equationResidual(end));

N = size(node,1);
udofNode = [node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2];
speedDof = sqrt(soln.u.^2+soln.v.^2);
speedNode = speedDof(1:N);

figs = gobjects(3,1);

figs(1) = figure('Visible',visible);
[~,bdFace] = findboundary3(elem);
trisurf(bdFace,node(:,1),node(:,2),node(:,3),speedNode,...
    'EdgeColor',[0.65 0.65 0.65],'FaceColor','interp');
axis equal tight;
view(35,24);
colorbar;
xlabel('x'); ylabel('y'); zlabel('z');
title('boundary speed: |U_h| on boundary');
exportfigure(figs(1),outputDir,'FirstOrder3Periodic-boundary-speed.png',exportPng);

ns = getoption(option,'sliceResolution',81);
[xs,ys] = meshgrid(linspace(0,4,ns),linspace(0,4,ns));
zs = 0.5*ones(size(xs));
[us,vs,speeds] = interpolatevelocity(udofNode,soln,xs,ys,zs);

figs(2) = figure('Visible',visible);
surf(xs,ys,zs,speeds,'EdgeColor','none');
axis equal tight;
view(2);
colorbar;
xlabel('x'); ylabel('y');
title('z-slice speed: |U_h| on z = 0.5');
hold on;
stride = max(1,round(ns/16));
quiver3(xs(1:stride:end,1:stride:end),...
        ys(1:stride:end,1:stride:end),...
        zs(1:stride:end,1:stride:end),...
        us(1:stride:end,1:stride:end),...
        vs(1:stride:end,1:stride:end),...
        zeros(size(us(1:stride:end,1:stride:end))),...
        'k');
hold off;
exportfigure(figs(2),outputDir,'FirstOrder3Periodic-mid-slice-speed.png',exportPng);

bedNode = abs(node(:,3)) < 10*eps;
bedElem = elem(all(bedNode(elem),2),:);
if isempty(bedElem)
    bedTri = [];
else
    bedTri = bedElem(:,[1 2 3]);
end

figs(3) = figure('Visible',visible);
if isempty(bedTri)
    scatter3(node(bedNode,1),node(bedNode,2),node(bedNode,3),36,...
        speedNode(bedNode),'filled');
else
    trisurf(bedTri,node(:,1),node(:,2),node(:,3),speedNode,...
        'EdgeColor',[0.65 0.65 0.65],'FaceColor','interp');
end
axis equal tight;
view(2);
colorbar;
xlabel('x'); ylabel('y');
title('basal speed: |U_h| on bed z = 0');
exportfigure(figs(3),outputDir,'FirstOrder3Periodic-bed-speed.png',exportPng);

if nargout == 0
    clear figs
end

end

function [us,vs,speeds] = interpolatevelocity(udofNode,soln,xs,ys,zs)
Fu = scatteredInterpolant(udofNode(:,1),udofNode(:,2),udofNode(:,3),...
    soln.u,'natural','none');
Fv = scatteredInterpolant(udofNode(:,1),udofNode(:,2),udofNode(:,3),...
    soln.v,'natural','none');
us = Fu(xs,ys,zs);
vs = Fv(xs,ys,zs);
speeds = sqrt(us.^2+vs.^2);
end

function exportfigure(fig,outputDir,fileName,exportPng)
if ~exportPng, return; end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
exportgraphics(fig,fullfile(outputDir,fileName),'Resolution',200);
end

function value = getoption(option,name,defaultValue)
if isfield(option,name)
    value = option.(name);
else
    value = defaultValue;
end
end
