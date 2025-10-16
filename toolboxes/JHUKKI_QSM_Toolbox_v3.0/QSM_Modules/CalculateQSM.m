%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09

%% Get variables
Params      = handles.Params;
GREMag      = handles.GREMag;
maskErode   = handles.maskErode;
freqMap     = handles.freqMap;

if prod(Params.sizeVol) > 1e7 % in case of ultra high resolution
    datatype = 'single';
else
    datatype = 'double';
end

% estimate noise level here
noise_level = calfieldnoise(GREMag(:,:,:,:,1), handles.GREPhaseRaw1, maskErode);

% Calculate base values
deltaB          = freqMap./(Params.gamma*Params.B0)*1e6;   % in ppm, field map
deltaB          = deltaB.*maskErode;                       
clear freqMap GREPhase GREPhaseRaw GREPhaseRaw1

deltaB = cast(deltaB, datatype);
GREMag = cast(GREMag, datatype);

%% switch for getting R2* and AutoRef (SFCR+0)
if Params.R2starFlag == 1
    % fitting R2*
    R2starFile = [Params.FileBaseName, '_R2star'];         
    if ~exist(fullfile(cd, [R2starFile, '.mat']), 'file') && size(GREMag,4)>2         
        fname1 = [Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(1))];
        BrainMaskFilename = [fname1, '_brain_mask.nii.gz'];
        
        nii = load_untouch_nii(BrainMaskFilename);      
        maskBET = permute(nii.img, [2,1,3]);            
        maskBET = maskBET > 0;
        
        % ------------- calculate R2* Maps and save file
        disp('Fitting R2* ...')
        tic
            R2starmap = R2star_ARLO(GREMag, Params.TEs, maskBET);        % ARLO faster            
        toc
        disp('Done.')
        
        R2starMax = 500;
        R2starmap(R2starmap < 0) = 0;
        R2starmap(R2starmap > R2starMax) = R2starMax;
        
        save(fullfile(cd, [R2starFile, '.mat']), 'R2starmap', 'Params', 'maskBET');
        saveNII(R2starmap, R2starFile, Params);

        clear maskBET nii fname1 R2starmap S0map R2starFile
        delete(gcp('nocreate'))    
    else
        disp('R2* fitting has been done already.')
    end
    
    if (Params.AutoRefFlag == 1) && (Params.B0 == 3)    % if AutoRef and 3T
       
        fileName = [ Params.FileBaseName, '_GREMag1'];
        if ~exist(fullfile(cd, [fileName, '.nii']), 'file')
            saveNII(GREMag(:,:,:,1).*1, fileName, Params, 1);
        end

        BrainMaskFilename = [fileName, '_brain_mask.nii.gz'];
        if ~exist(fullfile(cd, BrainMaskFilename), 'file')
            inputstring1 = [Params.FSLFolder, 'bet ', fileName, '.nii ' fileName, '_brain', ' -f 0.5 -g 0 -m' ];
            disp(inputstring1)
            system(inputstring1);
        end
        disp('bet done.')

        % do fast on GREMag 1 to improve CSF mask
        BETFilename = [fileName, '_brain'];   
        FASTFilename = [fileName, '_brain_seg.nii.gz'];
        if ~exist(fullfile(cd, FASTFilename), 'file') 
            inputstring1 = [Params.FSLFolder, 'fast -t 1 -n 2 -H 0.1 -I 4 -l 20.0 -o ', BETFilename, ' ', BETFilename];
            disp(inputstring1)
            system(inputstring1);
        end
        disp('fast done.')
        clear BETFilename FASTFilename BrainMaskFilename
    end    
 
end

if Params.AutoRefFlag == 1
    % -------------  load in R2* and fsl fast segmetnation
    R2starFile = [Params.FileBaseName, '_R2star'];     
    if (exist(fullfile(cd, [R2starFile, '.mat']), 'file') == 2)
        disp('R2* map exist.')
        temp = load(R2starFile);
        R2starMap = temp.R2starmap;
        clear temp
        R2starMap = cast(R2starMap, datatype);
    end

    GREMagSegFileFlag = 0;
    GREMagSegFile = [Params.FileBaseName, '_brain_mixeltype.nii.gz'];             % old version
    GREMagSegFile2 = [Params.FileBaseName, '_GREMag1_brain_mixeltype.nii.gz'];    % new version
    if (exist(fullfile(cd, GREMagSegFile), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile;            % .nii.gz
        GREMagSegFileFlag = 1;
    elseif (exist(fullfile(cd, GREMagSegFile(1:end-3)), 'file') == 2)
        GREMagSegFileTarget = GREMagSegFile(1:end-3);   % .nii
        GREMagSegFileFlag = 1;
    elseif (exist(fullfile(cd, GREMagSegFile2), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile2;
        GREMagSegFileFlag = 1;    
    elseif (exist(fullfile(cd, GREMagSegFile2(1:end-3)), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile2(1:end-3);   % .nii
        GREMagSegFileFlag = 1;
    end

    if GREMagSegFileFlag == 1
        disp('GREMagSeg file exists.')
        nii = load_untouch_nii(GREMagSegFileTarget);
        GREMagSeg = permute(nii.img, [2,1,3]);   
    end
end

% ----------------------
GREMag = GREMag./noise_level;

%% make DPWeight if needed
switch Params.QSMSolverDict{Params.QSMSolver} 
    case {'MEDI', 'SFCR'}
        DPWeight = zeros([Params.sizeVol, 1, Params.nDynamics]);
        if isfield(handles, 'DPWeight')
            DPWeight = handles.DPWeight;        % loaded weight may not be normalized
        else
            for dynamic_ind = 1:Params.nDynamics
                DPWeight(:,:,:,1,dynamic_ind) = sqrt(sum(single(GREMag(:,:,:,Params.echoNums,dynamic_ind)).^2, 4));
                temp = DPWeight(:,:,:,1,dynamic_ind);
                DPWeight(:,:,:,1,dynamic_ind) = DPWeight(:,:,:,1,dynamic_ind)./mean(temp(maskErode>0));  % normalization
                DPWeight(:,:,:,1,dynamic_ind) = DPWeight(:,:,:,1,dynamic_ind).*maskErode; 
            end
        end
        clear GREMag  
    otherwise
        clear GREMag
end

%% determin output file name
switch Params.QSMSolverDict{Params.QSMSolver} 
    case 'iLSQR'
        StringApp1 = '_chi_iLSQR';  WaitBarMsgloading = 'iLSQR';
    case 'TKD'
        StringApp1 = '_chi_TKD';    WaitBarMsgloading = 'TKD';        
    case 'iTKD'
        StringApp1 = '_chi_iTKD';   WaitBarMsgloading = 'iTKD';
    case 'MEDI'
        StringApp1 = '_chi_MEDI';   WaitBarMsgloading = 'MEDI';        
    case 'SFCR'
        if Params.AutoRefFlag ==1
            StringApp1 = '_chi_SFCR+0';   WaitBarMsgloading = 'SFCR+0';  
        else
            StringApp1 = '_chi_SFCR';   WaitBarMsgloading = 'SFCR';  
        end
    otherwise 
        error('unknown selection of QSM method.')
end

% naming convention
if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
    StringApp2 = '_NLSlope';
elseif Params.EchoAvg > 0
    StringApp2 = '_Avg';
else
    StringApp2 = '_Slope';
end

if (length(Params.echoNums) > 1)            
    StringApp4 = ['-', num2str(Params.echoNums(end))];
else
    StringApp4 = '';
end

StringApp3 = ['_echo', num2str(Params.echoNums(1))];

if Params.echoStep > 1      % For Single Echo, Params.echoStep=0;
    StringApp5 = ['_s', num2str(Params.echoStep)];
else
    StringApp5 = '';
end

outputFile = [Params.FileBaseName, StringApp1, StringApp2];     % version 1
% outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, StringApp4, StringApp5];

%% QSM solver
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && ...
        Params.saveOutput && ...
        ~Params.BgWasRemoved) % Re-DO QSM if new background removal was performed!!
    multiWaitbar( ['Loading previously ', WaitBarMsgloading,  ' calculated QSM...'], 'Busy');
    load(outputFile)
    multiWaitbar('CloseAll');
else            
    textWaitbar = ['Calculating QSM for '  num2str(Params.nDynamics), ' dynamics.'];
    multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
    chi_res = zeros(size(deltaB));
    
    switch Params.QSMSolverDict{Params.QSMSolver} 
        case 'iLSQR'        
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype); 
            
            for dynamic_ind = 1:Params.nDynamics
                
                % initial estimation using LSQR2
                disp('initial LSQR...')
                tol = 0.01;            
                chi_res0 = delta2chi_lsqr2(deltaB(:,:,:,1,dynamic_ind), D, tol, maskErode);                   
                chi_res0 = chi_res0.*maskErode;

                disp('estimating a-priori using iTKD ...')
                tol = 0.04;       
                [chi_ap] = delta2chi_iTKD(deltaB(:,:,:,1,dynamic_ind), D, maskErode, tol);    

                % further remove possible streaking artifact  if necessary            
                thresh_SAR = 0.15;
                [chi_ap, chi_SA_ap] = delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_ap, maskErode, thresh_SAR, chi_ap);  

                % Final estimate streaking artifact in LSQR2
                disp('final SAR ...')
                thresh_SAR = 0.15;    
                [chi_res(:,:,:,1,dynamic_ind), chi_SA] = delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_res0, maskErode, thresh_SAR, chi_ap);
                
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end

        case 'TKD' 
            thresh_tkd = 0.2;
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            
            for dynamic_ind = 1:Params.nDynamics
                chi_res(:,:,:,1,dynamic_ind) = delta2chi_tso(deltaB(:,:,:,1,dynamic_ind), D, thresh_tkd);
                chi_res(:,:,:,1,dynamic_ind) = chi_res(:,:,:,1,dynamic_ind).*maskErode;
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end

        case 'iTKD'    
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);            

            tol = 0.04;       
            for dynamic_ind = 1:Params.nDynamics
                chi_res0 = delta2chi_iTKD(deltaB(:,:,:,1,dynamic_ind), D, maskErode, tol);   
                thresh_SAR = 0.15;
                [chi_res(:,:,:,1,dynamic_ind), chi_SA] = ...
                    delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_res0, maskErode, thresh_SAR, chi_res0);      
                
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end
            
        case 'MEDI'     
            lambda = 2500;       
            edgePer = 3*0.3;    % Edge voxel percentage
            merit = 1;          % fine tuning
            
            for dynamic_ind = 1:Params.nDynamics
                
                [chi_res(:,:,:,1,dynamic_ind), regv, datav] = ...
                    delta2chi_MEDI(deltaB(:,:,:,1,dynamic_ind), Params, DPWeight(:,:,:,1,dynamic_ind), maskErode, lambda, merit, edgePer);  
                
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;               
            end
            
        case {'SFCR'}          %% modifed SFCR
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            D = ifftshift(D);
            
            % Parameters set for modified SFCR code (2&3)
            MagEdgePer = 0.3;  % Edge voxel percentage for magnitude 
            chiEdgePer = 0.35;

            lambdaSet.beta = 10;
            lambdaSet.lambda1_M = 500;      
            lambdaSet.lambda2_M = 0;
            lambdaSet.lambda1_S = 1000;     
            lambdaSet.lambda2_S = 0;
            
            thre_tkd = 0.18;                % 
            
            if (Params.AutoRefFlag == 1) && (exist('R2starMap', 'var') == 1)                
                lambdaSet.R2sThresh = 5;        % 5 Hz for extracting central CSF region for automatic CSF referencing               
                [CSFmask1] = CSFmaskThresh(R2starMap, lambdaSet.R2sThresh, maskErode, Params.voxSize);
                
                if (exist('GREMagSeg', 'var') == 1)
                    CSFmask2 = (GREMagSeg == 0) & maskErode;
                    disp('updating CSF mask based on R2* with FSL Segmentation')
                else
                    CSFmask2 = maskErode;
                end
                
                lambdaSet.maskSS = CSFmask1 & CSFmask2;
                lambdaSet.lambda2_M = lambdaSet.lambda1_M./5;
                lambdaSet.lambda2_S = lambdaSet.lambda1_S./5;
            end

            for dynamic_ind = 1:Params.nDynamics
                if strcmp(Params.QSMSolverDict{Params.QSMSolver}, 'SFCR')
                    [chi_res(:,:,:,1,dynamic_ind), SB_residual, SB_regularization] = ...
                        delta2chi_SFCR(deltaB(:,:,:,1,dynamic_ind), D, maskErode, DPWeight(:,:,:,1,dynamic_ind), ...
                        MagEdgePer, chiEdgePer, lambdaSet, thre_tkd, Params.AutoRefFlag);
                end
                                
                hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                HandleStopReconstruction;
            end
        otherwise
            error('unknown selection of QSM method.')
    end     
end         

%% save
switch Params.QSMSolverDict{Params.QSMSolver} 
    case 'iLSQR'
        save([outputFile '.mat'], 'chi_res', 'chi_SA', 'chi_ap', 'tol', 'thresh_SAR', 'maskErode');  
    case 'TKD'
        save([outputFile '.mat'], 'chi_res', 'thresh_tkd', 'maskErode');
    case 'iTKD'
        save([outputFile '.mat'], 'chi_res', 'chi_SA', 'chi_res0', 'thresh_SAR', 'tol', 'maskErode');
    case 'MEDI'
        save([outputFile '.mat'], 'chi_res', 'maskErode', 'Params', 'lambda', 'merit');  
    case 'SFCR'
        save([outputFile '.mat'], 'chi_res', 'maskErode', 'Params', 'MagEdgePer', ...
                'chiEdgePer', 'lambdaSet', 'SB_residual', 'SB_regularization');
    otherwise
        error('unknown selection of QSM method.')
end

multiWaitbar('CloseAll');
saveNII(squeeze(chi_res), outputFile, Params, 1);

% Save
handles.chi_res        = chi_res;
guidata(hObject, handles);

% Save Constants file
Params          = handles.Params;
save(handles.Params.QSMSettingsFile, 'Params');

% Display
handles = LoadImage(hObject, handles, handles.chi_res, 'Susceptibility (ppm)');

% Update panel
set(handles.PanelStep4, 'HighlightColor', [0 0.5 0]);
set(handles.ButtonShowQSM, 'Enable', 'On')

% Update table
UpdateTable(handles, 'Completed 4 of 4 - Finished!');

% Update buttons
set([handles.ButtonAddDataset handles.ButtonStartDatasets handles.ButtonLoadDataList], 'Enable', 'On')

% Enable buttons
set([handles.ButtonEditEchoes handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver handles.VarFSLThres], 'Enable', 'On');