%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09

%% Perform Phase Unwrapping
% Remove open waitbars!
wh=findall(0,'tag','TMWWaitbar');
delete(wh);

Params = handles.Params;
% ------------------------------------------------------------
% Disable buttons
set([handles.ButtonEditEchoes handles.VarFSLThres handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver], 'Enable', 'Off');

% Update panel
set([handles.PanelStep1 handles.PanelStep2 handles.PanelStep3 handles.PanelStep4], 'HighlightColor', [0 0 0]); % Black
set(handles.PanelStart, 'HighlightColor', [0 0.5 0]); % Green

% Get variables
GREPhase    = handles.GREPhase;

% Output
StringApp1 = ['_PhaseUnwrapped_echo' num2str(Params.echoNums(1))]; 

if(length(Params.echoNums) > 1)
    if Params.echoStep > 1
        StringApp2  = ['-' num2str(Params.echoNums(end)), '_s', num2str(Params.echoStep)];        
    else
        StringApp2  = ['-' num2str(Params.echoNums(end))];
    end
else
    StringApp2  = '';
end

switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod}
    case 'Path'
        StringApp3 = 'Path';        
    case 'Laplacian'
        StringApp3 = 'Lap';
    case 'NonlinearFit + Path'
        StringApp3 = 'NLFpath';
    otherwise
        error('Unknown unwrapping method.')
end

if Params.nDynamics == 1
    stringApp4 = '';    % default with just 1 dynamic
else
    stringApp4 = ['_dyn', num2str(Params.nDynamics)]; % 5D
end

outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, stringApp4];

% Switch
cd(Params.PathName);

% Check if exists
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && Params.saveOutput)
    multiWaitbar('Loading previously unwrapped phase...', 'Busy');
    load(outputFile)      
    multiWaitbar('CloseAll');
else
    textWaitbar = ['Unwrapping ' num2str(length(Params.TEs)) ' echoes, ', num2str(Params.nDynamics), ' dynamics'];
    multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
    
    switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
        case 'Laplacian'
            for dynamic_ind = 1:Params.nDynamics
                for echo_ind = 1:length(Params.TEs)
                    % Unwrapping all the echoes
                    if Params.phase2DprocFlag == 0  
                        % default is 3D
                        GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params);
                    else
                        % 2D phase unwrapping slice by slice
                        textWaitbar2D = 'Performing 2D phase unwrapping first';
                        multiWaitbar(textWaitbar2D, 0, 'Color', 'b' );

                        for sliceii = 1:Params.sizeVol(3)
                            GREPhase(:,:,sliceii,echo_ind,dynamic_ind) = phase_unwrap_laplacian_2D(GREPhase(:,:,sliceii,echo_ind,dynamic_ind), 0);  % no refVox
                            multiWaitbar(textWaitbar2D, sliceii/Params.sizeVol(3), 'Color', 'b' );
                        end            
                        GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params);
                    end
                    % Waitbar        
                    hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end                

        case 'NonlinearFit + Path'
            if length(Params.echoNums) > 2  % use nonlinear formula  
                GREMag = handles.GREMag;    % complex fitting need GREMag too   
                for dynamic_ind = 1:Params.nDynamics                     
                    [GREPhase(:,:,:,1,dynamic_ind), N_std, relres, ~] = Fit_ppm_complex(GREMag(:,:,:,Params.echoNums,dynamic_ind).*cos(GREPhase(:,:,:,Params.echoNums,dynamic_ind)) ...
                                                        + 1i*GREMag(:,:,:,Params.echoNums,dynamic_ind).*sin(GREPhase(:,:,:,Params.echoNums,dynamic_ind)));                                            
                    multiWaitbar( textWaitbar, 0.6*(dynamic_ind/Params.nDynamics));                                
                    GREPhase(:,:,:,1,dynamic_ind) = phase_unwrap_path_mex(double(GREPhase(:,:,:,1,dynamic_ind)));    % path based mex version, still in radian
                    hasCanceled = multiWaitbar(textWaitbar, dynamic_ind/Params.nDynamics);  

                    handles.DPWeight = 1./N_std;
                    handles.DPWeight(isinf(handles.DPWeight)) = 0;
                end
            else  % 1 or 2 echoes, otherwise just do pathbased unwrapping, change GREPhase to only one volume
                for dynamic_ind = 1:Params.nDynamics 
                    GREPhase(:,:,:,:,dynamic_ind) = double(GREPhase(:,:,:,:,dynamic_ind));
                    for echo_ind = 1:length(Params.echoNums)
                        GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind) = phase_unwrap_path_mex(GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind));
                        hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.echoNums))*(dynamic_ind/Params.nDynamics) );
                    end
                    if length(Params.echoNums) > 1
                        GREPhase(:,:,:,1,dynamic_ind) = find_slopeintercept_phasevstime_fast(GREPhase(:,:,:,Params.echoNums,dynamic_ind), [0,1]);
                    else
                        GREPhase(:,:,:,1,dynamic_ind) = GREPhase(:,:,:,Params.echoNums,dynamic_ind);
                    end
                end
            end
            GREPhase = GREPhase(:,:,:,1,:); % echo dimention collapse to 1
                        
        case 'Path'
            GREPhase = double(GREPhase);
            N = size(GREPhase);
            RefVox = [floor(N(1)/2), floor(N(2)/2), floor(N(3)/2)]; 

            for dynamic_ind = 1:Params.nDynamics
                for echo_ind = 1:length(Params.TEs)
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_path_mex(GREPhase(:,:,:,echo_ind,dynamic_ind));
                    hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics));
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = GREPhase(:,:,:,echo_ind,dynamic_ind) -  GREPhase(RefVox(1), RefVox(2), RefVox(3), echo_ind,dynamic_ind);
                end
            end

        otherwise
             error('Unknown unwrapping method.')

    end
    save([outputFile '.mat'], 'GREPhase', 'Params');
    multiWaitbar('CloseAll');
end

% Save
handles.GREPhase        = GREPhase;
guidata(hObject, handles);

% Display
handles = LoadImage(hObject, handles, handles.GREPhase, 'Phase (rad)');

% Update panel
set(handles.PanelStep1, 'HighlightColor', [0 0.5 0]);
set(handles.ButtonShowPhase, 'Enable', 'On')
% Update table
UpdateTable(handles, 'Completed 1 of 4');
