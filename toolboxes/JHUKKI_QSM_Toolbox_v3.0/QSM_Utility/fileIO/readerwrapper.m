function [GREMag, GREPhase, Params, handles] = readerwrapper(PathName, FileName, handles)
% function [GREMag, GREPhase, Params] = readerwrapper(PathName, FileName, handles)
%% Authors: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu
%
% readfilewrapper.m is the wrapper of readers for different file formats 
% created for QSM_Toolbox_v3
% updated 2019-07-09

[~,FileBaseName,FileExt] = fileparts(FileName);

if(strcmpi(FileExt,'.par'))
    %% If GRE data in Rec Format
    GREdataAll  = readrec_V4_4([PathName FileName], 1);  

    % Parameters
    Params      = ReadParFile([PathName FileName], handles.Params);     % extract basic scan parameters from .par file
    Params      = permuteParams(Params);                                

    GREMag = (GREdataAll(:,:,:,:,:,1,1,1,1));     % data type Mag, ncol, nrow, nslice,nechoes,ndynamics
    GREPhase = (GREdataAll(:,:,:,:,:,1,1,1,2));   % data type Phase same 5D array

    %% Slice orientation

    switch Params.sliceOri
        case 1
            %% TRANSVERSAL (Normal)
            % Do nothing
        case 3
            %% CORONAL
            % Not fully tested yet
        case 2
            % Not fully tested yet
        otherwise
            disp('error info in slice orientation')
    end

elseif(sum(strcmpi(FileExt,{'.DIC';'.IMA';'.DICOM'; '.dcm'; '.1'})) > 0)
    % check whether conventional/enhanced dicom
    disp('reading dicom info ...')
    dicomheader = dicominfo([PathName FileName]);
    disp('Done.')

    if ~isfield (dicomheader, 'PerFrameFunctionalGroupsSequence')
        %% conventional DICOM in PathName
        [GREPhase, GREMag, Params] = read_DICOM(PathName, handles.Params);
        GREMag = permute(GREMag, [2,1,3:length(size(GREMag))]);
        GREPhase = permute(GREPhase, [2,1,3:length(size(GREPhase))]);
        Params = permuteParams(Params);

    else
        %% enhanced DICOM
        disp('reading enhanced DICOM data, please wait ...')
        [GREdataAll, dicomheader] = dicomeread([PathName FileName]);
        disp('Done.')
        Params = readparamsfromdicom(dicomheader, handles.Params);    
        [Params.PathName, Params.FileBaseName, ~] = fileparts([PathName, FileName]);

        ndimsGREdataAll = ndims(GREdataAll);
        if ndimsGREdataAll == 5
            GREdataAll = permute(GREdataAll ,[1,2,3,5,4]);  % 3D x imagetype x echoes
            GREMag = GREdataAll(:,:,:,:,1);     % if multi-echo, Magnitude
            GREPhase = GREdataAll(:,:,:,:,2);   % Phase
        elseif ndimsGREdataAll == 4
            GREMag = GREdataAll(:,:,:,1);       
            GREPhase = GREdataAll(:,:,:,2);
        end

    end
    set(handles.VarB0,'String', Params.B0);    

    saveDICOM2mat = 1;
    if saveDICOM2mat == 1
        save([Params.FileBaseName, '.mat'], 'GREMag', 'GREPhase', 'Params');
    end


elseif(strcmpi(FileBaseName, 'method'))
    %% Bruker files
    %  !!!!!!!!!   not updated for a while, use with caution!!!!
    % Load params
    Params = ReadBrukerParams(PathName, 1, handles.Params);

    % Update labels
    set(handles.TextFileName, 'String', Params.FileBaseName);
    set(handles.VarB0, 'String', '11.7'); % GUESS

    % Load data
    %% Combine phase array data from each coil

    coilcombine_flag = 2;
    if coilcombine_flag == 1
        GREMag = read_2dseq_V1(PathName, 1);
        GREPhase = read_2dseq_V1(PathName, 2);
        GREMag = single(GREMag);
        GREPhase = single(GREPhase);

        % Coil combination PHASE
        GREPhase = mcpcvcr(GREMag, GREPhase, Params);        % Multi-channel phase combination, with mcpvcr

        % Coil combination MAGNITUDE
        if Params.nEchoes > 1   % multiecho data
            GREMag = sqrt(sum(GREMag.^2, 5));                         
            Params.sizeRecAll = Params.sizeRecAll(1:4);
        else
            GREMag = sqrt(sum(GREMag.^2, 4));                  
            Params.sizeRecAll = Params.sizeRecAll(1:3);
        end            
    else
        GREMag = read_2dseq_V1(PathName, 1);
        GREPhase = read_2dseq_V1(PathName, 2);
        GREMag = single(GREMag);
        GREPhase = single(GREPhase);
    end

elseif (strcmpi(FileExt,'.mat'))
    % matlab .mat file, need to have both GREMag, GREPhase, Params
    S = load([PathName FileName]);

    if isfield(S.Params, 'B0')
        set(handles.VarB0,'String', S.Params.B0);     
    end

    GREMag = S.GREMag;
    GREPhase = S.GREPhase;
    Params = handles.Params;

    Params.sizeVol = S.Params.sizeVol;
    Params.fov = S.Params.fov;
    Params.voxSize = S.Params.voxSize;

    Params.B0 = S.Params.B0;
    Params.TR = S.Params.TR;
    Params.nEchoes = S.Params.nEchoes;
    Params.TEs = S.Params.TEs;

    Params.TAng = S.Params.TAng;

    Params.PathName = PathName;
    Params.FileBaseName = FileBaseName;

    if isfield(S.Params, 'nDynamics')
        Params.nDynamics = S.Params.nDynamics;
    else
        Params.nDyanmics =  1;
    end

    if isfield(S.Params, 'sliceOri')
        Params.sliceOri = S.Params.sliceOri;
        Params.ang = S.Params.ang;
        Params.AngAP = S.Params.AngAP;
        Params.AngFH = S.Params.AngFH;
        Params.AngRL = S.Params.AngRL;            
    end

    if isfield(S.Params, 'Tsom')
        Params.TAnginv = S.Params.TAnginv;

        Params.Tpom = S.Params.Tpom;
        Params.Tpominv = S.Params.Tpominv;

        Params.Tsom = S.Params.Tsom;
        Params.Tsominv = S.Params.Tsominv;
    end

else
    error('Sorry, we can not open this type of file...');
end