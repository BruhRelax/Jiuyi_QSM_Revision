%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated 2019-07-09 X.L.

%% Start processing of multiple datasets
% Save Constants file
Params          = handles.Params;
save(handles.Params.QSMSettingsFile, 'Params');

% Get table data
tableData = get(handles.TableDatasets, 'Data');

% Disable buttons
set([handles.ButtonEditEchoes handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver], 'Enable', 'Off');

% Go through all the data
for c = 1:size(tableData, 1)
    if(strcmpi(tableData{c,2}, handles.textReadyLoad) || strcmpi(tableData{c,2}, handles.textReconstructing))
        % Update number and text
        handles.CurrentDataset = c;
        UpdateTable(handles, 'Reconstructing...');
        
        % Filename
        FileFullName = tableData{c,1};
        
        % Extract info
        [PathName, FileBaseName, FileExt] = fileparts(FileFullName);
        
        % Switch to that folder
        cd(PathName);
        
        % readin data
        [GREMag, GREPhase, Params, handles] = readerwrapper([PathName, filesep], [FileBaseName, FileExt], handles);
        
        % Save to workspace
        handles.GREPhase    = single(GREPhase);       
        handles.GREPhaseRaw = single(GREPhase);              % for backup raw phase
        handles.GREPhaseRaw1 = single(GREPhase(:,:,:,1,1));    % just for display using the first echo/dynamic, 3D
        handles.GREMag      = single(GREMag);
        handles.slice2disp  = floor(Params.sizeVol(3)/2);
        handles.Params      = Params;
        
        % clear redundant
        clear GREMag GREPhase

        handles = LoadImage(hObject, handles, handles.GREPhaseRaw1, 'Phase (rad)');
    
        % Update info
        ShowImageInfo;
        
        % Enable slider
        set(handles.ImageSlider,'Enable','On')
        set(handles.ImageSlider,'Value',handles.slice2disp)
        set(handles.ImageSlider,'Min',1)
        set(handles.ImageSlider,'Max',handles.Params.sizeVol(3))
        
        %% LETS DO THISSSSS
        PerformUnwrapping;
        CreateBrainMask;
        RemoveBackground;
        CalculateQSM;
    end
    
end

