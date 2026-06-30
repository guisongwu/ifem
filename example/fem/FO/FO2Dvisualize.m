function fig = FO2Dvisualize(option,soln,eqn,info,node,elem,bdFlag)
%% FO2DVISUALIZE visualize the 2D x-z FO solution.
%
% Example:
%
%   FO2Dvisualize
%
% Reuse an existing solve:
%
%   [soln,eqn,info,node,elem,bdFlag] = FirstOrderP2(option);
%   FO2Dvisualize(option,soln,eqn,info,node,elem,bdFlag)

if nargin < 1, option = struct; end
visible = getoption(option,'visible','on');
exportPng = getoption(option,'exportPng',false);
outputDir = getoption(option,'outputDir','fig');

if nargin < 2 || isempty(soln)
    solveOption = option;
    solveOption.printlevel = getoption(option,'printlevel',1);
    [soln,eqn,info,node,elem,bdFlag] = FirstOrderP2(solveOption);
elseif nargin < 7
    error('iFEM:FO2DVisualizationInput',...
        'Provide soln, eqn, info, node, elem, and bdFlag together.');
end

fprintf('FO2D visualization: dof=%d, it=%d, residual=%.3e\n',...
    length(soln.u),info.itStep,info.equationResidual(end));

N = size(node,1);
dofNode = [node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2];
uNode = soln.u(1:N);

fig = figure('Visible',visible);
trisurf(elem,node(:,1),node(:,2),uNode,...
    'EdgeColor',[0.7 0.7 0.7],'FaceColor','interp');
view(2);
axis equal tight;
colorbar;
xlabel('x');
ylabel('z');
title('2D FO x-z cross-section: u(x,z)');

hold on;
nsx = getoption(option,'vectorNx',17);
nsz = getoption(option,'vectorNz',7);
slope = getgeometryfield(info,option,'slope',0);
H = getgeometryfield(info,option,'H',max(node(:,2))-min(node(:,2)));
xMin = min(node(:,1));
xMax = max(node(:,1));
[xq,qq] = meshgrid(linspace(xMin,xMax,nsx),linspace(0,H,nsz));
zq = qq-slope*xq;
Fu = scatteredInterpolant(dofNode(:,1),dofNode(:,2),soln.u,'natural','none');
uq = Fu(xq,zq);
quiver(xq,zq,uq,zeros(size(uq)),'k');
hold off;

if exportPng
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    exportgraphics(fig,fullfile(outputDir,...
        'FO2D-cross-section-u.png'),'Resolution',200);
end

if nargout == 0
    clear fig
end

end

function value = getoption(option,name,defaultValue)
if isfield(option,name)
    value = option.(name);
else
    value = defaultValue;
end
end

function value = getgeometryfield(info,option,name,defaultValue)
if isfield(info,'geometry') && isfield(info.geometry,name)
    value = info.geometry.(name);
elseif isfield(option,name)
    value = option.(name);
else
    value = defaultValue;
end
end
