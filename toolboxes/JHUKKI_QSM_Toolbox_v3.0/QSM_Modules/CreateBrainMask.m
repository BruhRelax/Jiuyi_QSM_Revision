%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09

%% Get variables
Params      = handles.Params;
GREPhase    = handles.GREPhase;
GREMag      = handles.GREMag;
GREPhaseRaw = handles.GREPhaseRaw;

% Output
if(length(Params.SaveEcho) > 1)
    % if using multiple echoes for masking
    outputFile  = [Params.FileBaseName '_brain_mask_echo' num2str(Params.SaveEcho(1)) '-' num2str(Params.SaveEcho(2)) '_r' num2str(Params.ErodeRadius) '_t' num2str(Params.FSLThreshold)];
else
    outputFile  = [Params.FileBaseName '_brain_mask_echo' num2str(Params.SaveEcho) '_r' num2str(Params.ErodeRadius) '_t' num2str(Params.FSLThreshold)];
end


%% Check if exists
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && Params.saveOutput)
    multiWaitbar('Loading previously created mask...', 'Busy');
    load([outputFile, '.mat'])
    multiWaitbar( 'CloseAll' );
    % Mark that we loaded mask
    handles.Params.MaskWasMade = 0;
else
    %% DOES FSL BET EXISITS?
    if(exist(Params.FSLFolder, 'dir') == 7)
        %% New mask
        textWaitbar = 'Creating brain mask...';
        multiWaitbar( textWaitbar, 0 );
        
        % First save the magnitude data
        for i = 1:length(Params.SaveEcho)
            fileName = [ Params.FileBaseName, '_GREMag', num2str(Params.SaveEcho(i))];
            
            % Remove previous files of FSL output OTHERWSIE FSL WILL NOT RUN
            delete([fileName '_brain.*']);
            delete([fileName '_brain_mask.*']);
            
            % Save new ones            
            saveNII(GREMag(:,:,:,Params.SaveEcho(i), 1).*1, fileName, Params, 1);  
        end

        %% Using BET to extract BrainMask
        fname1 = [ Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(1))];
        inputstring1 = [Params.FSLFolder, 'bet ', fname1, '.nii ' fname1, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];
        multiWaitbar( textWaitbar, 0.2 );
        system(inputstring1);
        
        % Save
        output = [fname1, '_brain_mask.nii.gz'];
        
        % Did it do anything?
        if(exist(fullfile(cd, output), 'file') ~= 2)
            errordlg('FSL BET did not run and no previous mask was found.. Reconstruction stopped');
            multiWaitbar( 'CloseAll' );
            return;
        end
        
        nii = load_untouch_nii(output);             
        maskBET1 = permute(nii.img, [2,1,3]);       
        maskBET1 = maskBET1 > 0;
        maskErode = maskBET1;
                
        if length(Params.SaveEcho) > 1
            % Second mask at later echo, needs to be more conservative
            % using BET
            fname2 = [ Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(2))];
            inputstring1 = [Params.FSLFolder, 'bet ', fname2, '.nii ' fname2, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];
            multiWaitbar( textWaitbar, 0.5 );
            system(inputstring1);
            
            % Combine BET masks
            multiWaitbar( textWaitbar, 0.7 );
            output = [fname2, '_brain_mask.nii.gz'];
            
            nii = load_untouch_nii(output);
            maskBET2 = permute(nii.img, [2,1,3]);
            maskBET2 = maskBET2 > 0;
             
            % V2
            zstart = 1;
            zend = 0.8*size(maskErode, 3);    
            ystart = 1; 
            yend1 = floor(size(maskErode, 1)*0.8); 

            xstart = floor(size(maskErode,2)*0.1);
            xend = floor(size(maskErode,2)*0.9);

            maskErode(ystart:yend1, xstart:xend,zstart:zend) = ...
               maskErode(ystart:yend1, xstart:xend, zstart:zend) & maskBET2(ystart:yend1, xstart:xend, zstart:zend);
        end

        % MaskOut Unreliable Phase
        if Params.nEchoes > 1    
            temp = GREPhaseRaw(:,:,:,:,1);  
            mask_intrinsic = (sum(abs(temp - min(temp(:))) < 10*eps, 4) == (Params.nEchoes-1)); % air/bone mask
        else
            mask_intrinsic = zeros(size(maskErode));
        end
        
        disp('estimating unreliable phase ...')
        TV = TVOP;

        switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod}
            case {'Path', 'NonlinearFit + Path'} 

                mask_unrelyPhase = zeros(size(maskErode));
                
                if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
                    echoNumPick = 1;  
                else                                    
                    echoNumPick = Params.echoNums;
                end
                
                for echoNumPickii = 1:length(echoNumPick)                          
                    tempGrad = TV*(GREPhase(:,:,:,echoNumPick(echoNumPickii),1).*maskErode);
                    mask_unrelyPhase_temp = sum(abs(tempGrad) >= pi, 4) > 2;   % residual wraps
                    mask_unrelyPhase_temp = imdilate(mask_unrelyPhase_temp, strel('disk', 1));
                    tempMask = ~mask_unrelyPhase_temp.*maskErode;              

                    CC = bwconncomp(tempMask, 6);       
                    numPixels = cellfun(@numel, CC.PixelIdxList);
                    [biggest, idx] = max(numPixels);
                    tempMask(CC.PixelIdxList{idx}) = 10;        
                    mask_unrelyPhase = mask_unrelyPhase | (tempMask ~= 10);
                end
         
            case 'Laplacian' % Laplacian based
                % if using iRSHARP can keep the original BET mask
                switch Params.BgRemovalMethodsDict{Params.BgRemoval}
                    case {'VSHARP','PDF','LBV+VSHARP'}
                        if Params.B0 == 3 && length(Params.TEs)>1
                            [~, echoNumPick] = min(abs(Params.TEs - 30e-3));       
                            unrelyPhaseThresh = 0.4;     % empirical number 

                        elseif Params.B0 == 7 && length(Params.TEs)>1
                            [~, echoNumPick] = min(abs(Params.TEs - 18e-3));       
                            unrelyPhaseThresh = 0.4;      
                        else
                            echoNumPick = 1;            
                            unrelyPhaseThresh = 0.4;    
                        end             
                        mask_unrelyPhase = create_mask_unrelyPhase(GREPhase(:,:,:,echoNumPick,1), unrelyPhaseThresh);       
                    case {'iRSHARP'}
                        mask_unrelyPhase = zeros(size(maskErode));
                    otherwise
                        error('unknown option.')
                end
            
            otherwise
                error('Unknown unwrapping method.');
        end
        
        mask_unrelyPhase = (mask_unrelyPhase | mask_intrinsic);        
        maskErode = (maskErode.*(maskErode - mask_unrelyPhase) > 0);
        maskErode = imerode3dslice(maskErode, strel('disk', 1));
        maskErode = imdilate3dslice(maskErode, strel('disk', 1)); 
        maskErode = imfill3(maskErode);     
        disp('Done')   
        
        multiWaitbar(textWaitbar, 0.8);        
        % Erosion
        if max(Params.voxSize) > 4*min(Params.voxSize)
            maskErode = imerode3dslice(maskErode, strel('disk', double(floor(Params.ErodeRadius./min(Params.voxSize)))));          
        else
            maskErode = imerode3(maskErode, floor(Params.ErodeRadius./min(Params.voxSize)), 1);          
        end
        maskErode = maskErode.*(maskErode - mask_intrinsic) > 0;
        
        % save maskBET
        multiWaitbar(textWaitbar, 0.9 );
        save([outputFile '.mat'], 'maskErode')
        
        % Save NIFTI
        saveNII(maskErode.*1, outputFile, Params, 1);    
        clear maskBET1 maskBET2 maskBET nii
        multiWaitbar( 'CloseAll' );

        % Mark that we made new mask and so we need to do new BG/QSM
        handles.Params.MaskWasMade = 1;
    else
        %% Backup in case FSL does not work!!
        % Is this correct threshold?
        isCorrect = false;
        while(~isCorrect)
            % Do the first masking
            [maskErode,~,~] = create_mask_erode(GREMag(:,:,:,Params.SaveEcho(1),1).*1, Params.ErodeRadius, Params.MaskThreshold);
            
            % Show
            handles.CurrentImage    = GREPhase(:,:,:,Params.SaveEcho(1),1).*maskErode;
            isNewImage = 1;
            ShowImage;
            drawnow;
            
            % Ask
            correctCheck = questdlg(sprintf('FSL BET not found.. Using back-up function.\nMASKING QUALITY IS BEST WHEN USING FSL BET!\nIs the brain correctly extracted or should the threshold be adapted?'),...
                'Masking','Correct','Change threshold','Stop reconstuction','Correct');
            
            % Switch
            switch correctCheck
                case 'Correct'
                    isCorrect = true;
                case 'Change threshold'
                    newThres = inputdlg('Enter new threshold (1-100):', 'Change threshold', [1,50], {num2str(Params.MaskThreshold)});
                    Params.MaskThreshold = str2double(newThres{:});
                case 'Stop reconstuction'
                    hasCanceled = true;
                    HandleStopReconstruction;
                otherwise
                    hasCanceled = true;
                    HandleStopReconstruction;
            end
        end
        % Mark that we made new mask and so we need to do new BG/QSM
        handles.Params.MaskWasMade = 1;
        % Not saving..
    end
    
    
end

% Save to workspace
handles.maskErode = maskErode;

% Show
if Params.EchoAvg > 0
    handles = LoadImage(hObject, handles, GREPhase(:,:,:,Params.echoNums(1),1).*maskErode, 'Phase (rad)');
else
    handles = LoadImage(hObject, handles, GREPhase(:,:,:,1,1).*maskErode, 'Phase (rad)');
end

% Save
guidata(hObject, handles);

% Update panel
set(handles.PanelStep2, 'HighlightColor', [0 0.5 0]);
set(handles.ButtonShowMask, 'Enable', 'On')
% Update table
UpdateTable(handles, 'Completed 2 of 4');