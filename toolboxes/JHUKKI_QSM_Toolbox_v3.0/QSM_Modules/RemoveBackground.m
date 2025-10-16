%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09

%% Get variables
Params      = handles.Params;
GREPhase    = handles.GREPhase;
GREMag      = handles.GREMag;
maskErode   = handles.maskErode;
GREPhaseRaw = handles.GREPhaseRaw;

if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
    StringApp_Avg = 'NLSlope';    
elseif Params.EchoAvg > 0
    StringApp_Avg = 'avg';
else
    StringApp_Avg = 'LSlope';
end

%% determined output file names
switch Params.BgRemovalMethodsDict{Params.BgRemoval}
    case 'VSHARP'
        %% Use V-SHARP algorithm, Average afterward
        StringApp0 = 'SMV';
        StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg];
        WaitBarMsgloading = ['V-SHARP ', StringApp_Avg];

    case 'PDF'
        %% DIPOLE FITTING
        StringApp0 = 'PDF';
        StringApp1 = ['_', StringApp0, '_', StringApp_Avg];
        WaitBarMsgloading = ['PDF ', StringApp_Avg];

    case 'LBV+VSHARP'
        %% using LBV+VSHARP
        StringApp0 = 'LBVSMV';
        StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg]; 
        WaitBarMsgloading = ['LBV+SHARP ', StringApp_Avg];

    case 'iRSHARP'
        %% iRSHARP 
        StringApp0 = 'iRSHARP';
        StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg]; 
        WaitBarMsgloading = ['iRSHARP ', StringApp_Avg];
        
    otherwise
        error('unknown option.')
end

%% naming conventions
if (length(Params.echoNums) > 1)            
    StringApp2 = '';
    StringApp4 = ['-', num2str(Params.echoNums(end))];
else
    StringApp2 = '';
    StringApp4 = '';
end

StringApp3 = ['_echo', num2str(Params.echoNums(1))];

if Params.echoStep > 1      % For Single Echo, Params.echoStep=0;
    StringApp5 = ['_s', num2str(Params.echoStep)];
else
    StringApp5 = '';
end

outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, StringApp4, StringApp5];
    
% normalized phase by 2*pi and TEs
if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
    % normalized slope
    if length(Params.echoNums) > 1
        GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums(2)) - Params.TEs(Params.echoNums(1))));  % in Hz
    else  
        GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums)));  % in Hz
    end
    Params.echoNums = 1;   
    
elseif (Params.EchoAvg > 0)
    for selectedEcho = 1:length(Params.TEs)   
        TEsSE = Params.TEs(selectedEcho);
        GREPhase(:,:,:,selectedEcho,:) = GREPhase(:,:,:,selectedEcho,:)./(2*pi*TEsSE);  % in Hz
    end 
else
    % if multi-echo, fit the linear slope
    normTEs = 0:(length(Params.echoNums)-1);
    for dynamic_ind = 1:Params.nDynamics
        GREPhase(:,:,:,1,dynamic_ind) = find_slopeintercept_phasevstime_fast(GREPhase(:,:,:,Params.echoNums,dynamic_ind), normTEs);
    end
    
    GREPhase = GREPhase(:,:,:,1,:); 
    GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums(2)) - Params.TEs(Params.echoNums(1))));
    Params.echoNums = 1;     
end

% Set selected echo for weighting
selectedEcho = Params.echoNums(1);

%% checking file existance
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && ...
        Params.saveOutput && ...
        ~Params.MaskWasMade) % Re-DO background removal for new mask
    multiWaitbar( ['Loading previously ', WaitBarMsgloading, ' calculated field...'], 'Busy');
    load([outputFile, '.mat']);
    multiWaitbar('CloseAll');
    % Mark that we loaded and no new QSM calculation is required
    handles.Params.BgWasRemoved = 0;
else
    %% Removing background gradient
    textWaitbar = ['Removing Background field for '  num2str(Params.nDynamics), ' dynamics.'];
    multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
    freqMap = zeros([Params.sizeVol, 1, Params.nDynamics]); 

    switch Params.BgRemovalMethodsDict{Params.BgRemoval}
        case 'VSHARP'
            %% Use V-SHARP algorithm
            clear GREMag
            for dynamic_ind = 1:Params.nDynamics
                if prod(Params.sizeVol)*floor(Params.SHARPradius./min(Params.voxSize)) > 1e8  % in case ultra-high res
                    % V-SHARP MULTI - image space version
                    [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = SHARP_Adaptive_Multi(GREPhase(:,:,:,:,dynamic_ind), maskErode, Params.SHARPradius, Params.thresh_tsvd, Params, handles);
                else            
                    % VSHARP k-space verion
                    radiusStep = min(Params.voxSize);
                    radiusArray = radiusStep:radiusStep:Params.SHARPradius;
                    [freqMap(:,:,:,1,dynamic_ind), dpfield_fit, mask_eval] = VSHARP_k(GREPhase(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles);
                end

                if Params.phase2DprocFlag == 1                      
                    if sum(mod(Params.sizeVol, 2)) > 0 
                        textWaitbar2D = 'Performing 2D V-SHARP';
                        multiWaitbar(textWaitbar2D, 0, 'Color', 'b' );

                        for sliceii = 1:Params.sizeVol(3)
                            [freqMap(:,:,sliceii,1,dynamic_ind), dpfield_fit2D] = SHARP_Adaptive_Multi2D(freqMap(:,:,sliceii,:,dynamic_ind), maskErode(:,:,sliceii,:), ...
                                                    Params.SHARPradius, Params.thresh_tsvd, Params);
                            multiWaitbar(textWaitbar2D, sliceii/Params.sizeVol(3), 'Color', 'b' );
                        end
                    else
                        [freqMap(:,:,:,1,dynamic_ind), ~, mask_eval] = VSHARP2D_k(freqMap(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles);                       
                    end
                end    
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end

            if exist('mask_eval', 'var')
                maskErode = mask_eval;  
            end

        case 'PDF'
            %% DIPOLE FITTING       
            for dynamic_ind = 1:Params.nDynamics
                if isfield(handles, 'DPWeight')
                    DPWeight = handles.DPWeight;
                else
                    DPWeight = single(GREMag(:,:,:,selectedEcho,dynamic_ind));
                end

                if Params.EchoAvg > 0
                    GREPhase(:,:,:,1,dynamic_ind) = mean(GREPhase(:,:,:,Params.echoNums,dynamic_ind), 4);
                 end

                [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = dipolefit_v2(GREPhase(:,:,:,1,dynamic_ind), maskErode, DPWeight, Params);

                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end

        case 'LBV+VSHARP'
            clear GREMag

            for dynamic_ind = 1:Params.nDynamics

                Params.lbv.tol = 0.01;      % default
                Params.lbv.depth = -1;      % default
                Params.lbv.peel = 0;        % similar to mask erosion

                % do LBV echo by echo
                for echoii = 1:length(Params.echoNums)
                    GREPhase(:,:,:,Params.echoNums(echoii),dynamic_ind) = LBV(GREPhase(:,:,:,Params.echoNums(echoii),dynamic_ind), maskErode, Params.sizeVol, Params.voxSize, ...
                                                            Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
                end

                % VSHARP k-space verion
                radiusStep = min(Params.voxSize);
                radiusArray = radiusStep:radiusStep:Params.SHARPradius;
                [freqMap(:,:,:,1,dynamic_ind), dpfield_fit, mask_eval] = VSHARP_k(GREPhase(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles); 

                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end
            maskErode = mask_eval;
            
        case 'iRSHARP'
            clear GREMag
            for dynamic_ind = 1:Params.nDynamics
                Params.iRSHARP_C = 0.25;    % add in Params to pass to iRHSARPv1 
                [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = iRSHARPv1(GREPhase(:,:,:,:,dynamic_ind), GREPhaseRaw(:,:,:,:,dynamic_ind), maskErode, Params, handles);
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end
            
        otherwise
            error('unknown option.')
    end

    % Save MAT file
    save([outputFile '.mat'], 'freqMap', 'Params', 'maskErode');     

    % Save NIFTI file
    saveNII(squeeze(freqMap), outputFile, Params, 1);    

    % Ready
    multiWaitbar('CloseAll');            
    % Mark that we made a new calculation and so we need a new QSM
    handles.Params.BgWasRemoved = 1;

end
    

%% Save and update display
handles.freqMap        = freqMap;
handles.maskErode      = maskErode;
guidata(hObject, handles);

% Display
handles = LoadImage(hObject, handles, handles.freqMap, 'Frequency (Hz)');

% Update panel
set(handles.PanelStep3, 'HighlightColor', [0 0.5 0]);
set(handles.ButtonShowBg, 'Enable', 'On')
% Update table
UpdateTable(handles, 'Completed 3 of 4');