%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-06-20

%% Show Image info in GUI
% Format
Params.voxSize = single(Params.voxSize);

% Update values in the information table
set(handles.TextDimensions, 'String', sprintf('Dimensions:  %d x %d x %d px', Params.sizeVol(1), Params.sizeVol(2), Params.sizeVol(3)))
set(handles.TextVoxelSize, 'String', sprintf('Voxel size:  %0.3g x %0.3g x %0.3g mm', Params.voxSize(1), Params.voxSize(2), Params.voxSize(3)))
set(handles.TextEchoes, 'String', sprintf('Echoes:  %d', Params.nEchoes))
set(handles.TextDynamics, 'String', sprintf('Dynamics:  %d', Params.nDynamics))

% How many echoes?
if(length(Params.echoNums) == 1)
    % One echoe
    set(handles.TextTEs, 'String', sprintf('TE:  %0.3g ms', Params.TEs(1)*1000))
else
    set(handles.TextTEs, 'String', sprintf('TE''s: %0.3g - %0.3g ms  (delta-TE: %0.3g ms)', Params.TEs(1)*1000, Params.TEs(end)*1000, (Params.TEs(2)-Params.TEs(1))*1000))    
end
