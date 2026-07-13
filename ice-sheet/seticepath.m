function seticepath()
%%SETICEPATH Put ice-sheet model directories first on the MATLAB path.
%
% Run iFEM's setpath once from the repository root, then run this helper
% before working in ice-sheet examples.  It keeps the renamed ice-sheet
% model directories ahead of same-named files elsewhere in iFEM while
% excluding documentation, outputs, and archived experiments.

root = fileparts(mfilename('fullpath'));
excluded = {'docs','doc','output','output_eps','archive','unused','copy'};

% Remove excluded ice-sheet paths that may already have been added by setpath.
currentPaths = strsplit(path,pathsep);
for k = 1:numel(currentPaths)
    if isunderroot(currentPaths{k},root) && hascomponent(currentPaths{k},excluded)
        rmpath(currentPaths{k});
    end
end

% Add only active source and test directories, with ice-sheet taking priority.
candidatePaths = strsplit(genpath(root),pathsep);
keep = ~cellfun('isempty',candidatePaths);
for k = find(keep)
    keep(k) = ~hascomponent(candidatePaths{k},excluded);
end
addpath(strjoin(candidatePaths(keep),pathsep),'-begin');

end

function tf = isunderroot(folder,root)
tf = strcmp(folder,root) || startsWith(folder,[root filesep]);
end

function tf = hascomponent(folder,names)
components = strsplit(folder,filesep);
tf = any(ismember(components,names));
end
