function [GREPhase, GREMag,Params] = read_DICOM(DICOMdir, Params, verbose)
% function [GREPhase, GREMag,Params] = read_DICOM(DICOMdir, Params, verbose)
%
%% Author: Xu Li
%% EDITED BY JIRI VAN BERGEN - 2014
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% read in DICOM format SWI images (Mag or Phase)
%
% input:
%        DICOMdir: DICOM directory
%        Params: The basic params from Constants.m
% output:
%        data: image data
%        Params.voxSize
%        Params.sizeVol
%        Params.B0: main field strength
%        Params.TE echo time (array if multi-echo, TE1 otherwise)
%        Params.TAng
%        Params.TR
%
% created on 2013-11-15, modified from DICOM reader from MEDI toolbox
% 
% Updated on 2014-10-10 for old Philips DICOM format, X.L.
% Updated on 2014-10-22 for Philips DICOM format, X.L.
% Updated on 2015-04-07 for SIEMENS uncombined data with ASCII code, X.L.
% Updated 2017-04-19, X.L. no round TEs
% Updated 2019-02-15, changed to original 5D/6D data format
% Updated 2019-04-22, added in support for GE data, X.L.
% Updated 2019-05-20, added in support for TOSHIBA/Canon data, X.L.

if nargin < 1
    % uigetdir get the DICOMdir
    DICOMdir = uigetdir(pwd, 'Select DICOM directory');
    if DICOMdir == 0
        error('DICOM folder not selected.')
    end
    verbose = 1;
elseif nargin < 3
    verbose = 1;
end

% Fancy waitbar
textWaitbar = 'Reading in DICOM-files';
multiWaitbar(textWaitbar, 0);

% Basics
cd(DICOMdir);
[Params.PathName, Params.FileBaseName,~] = fileparts(pwd);
Params.PathName = fullfile(Params.PathName, Params.FileBaseName);

% Get files
filelist = dir(DICOMdir);

ii=1;
while ii<=length(filelist)
    if filelist(ii).isdir==1 || ~isDICOM(filelist(ii).name)
        filelist = filelist([1:ii-1 ii+1:end]);   % skip folders
    else
        ii=ii+1;
    end
end

filenameBase = [DICOMdir '/' filelist(1).name];

info = dicominfo(filenameBase);

% check manufacture
if strfind(info.Manufacturer, 'SIEMENS')
    Manufacturer = 1;
elseif strfind(info.Manufacturer, 'Philips')
    Manufacturer = 2;
elseif strfind(info.Manufacturer, 'GE')
    Manufacturer = 3;
elseif strfind(info.Manufacturer, 'TOSHIBA')
    Manufacturer = 4;
else
    error('unknow Manufacturer ... ')
end

% basic information
Params.sizeVol(1) = single(info.Width);
Params.sizeVol(2) = single(info.Height);

Params.voxSize(1) = single(info.PixelSpacing(1));       % mm
Params.voxSize(2) = single(info.PixelSpacing(2));
Params.voxSize(3) = single(info.SliceThickness);        % mm

% CF = info.ImagingFrequency*1e6;                    % central frequency, Hz
Params.B0 = info.MagneticFieldStrength;                 % in T
Params.TR = info.RepetitionTime*1e-3;                   % TR in sec

% Angulation matrix
Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
Params.TAng = [Affine2D, cross(Affine2D(:,1), Affine2D(:,2))];

% Dynamics imaging
if isfield(info, 'NumberOfTemporalPositions')        % Philips/TOSHIBA
    Params.DynamicNum = info.NumberOfTemporalPositions;
else
    Params.DynamicNum = 1;
end

% coil element combination
if ~isfield(Params, 'coilNum')                              % Params not set yet
    Params.coilNum = 1;                                     % default, single coil or combined
    if isfield(info, 'Private_0051_100f') && ~isfield(Params, 'coilName')         % SIEMENS    
        Params.coilName(1) = cellstr(char((info.Private_0051_100f(:))'));        % coil name, could be ASCII code
    end
end

% Number of Echoes
minSlice = 1e10;
maxSlice = -1e10;
NumEcho = 0;
for ii = 1:length(filelist)
    info = dicominfo([DICOMdir '/' filelist(ii).name]);
    if info.SliceLocation<minSlice                  % find lowest slice
        minSlice = info.SliceLocation;              % in mm
        minLoc = info.ImagePositionPatient;
    end
    if info.SliceLocation>maxSlice                  % find highest slice
        maxSlice = info.SliceLocation;
        maxLoc = info.ImagePositionPatient;
    end
    if info.EchoNumber>NumEcho                      %
        NumEcho = info.EchoNumber;
    end
    
    if isfield(info, 'Private_0051_100f')              % get all the coil names
        if isempty(find(strcmpi(Params.coilName, cellstr(char((info.Private_0051_100f(:))'))) > 0, 1, 'first'))  % new coil name
            Params.coilName(length(Params.coilName)+1) = cellstr(char((info.Private_0051_100f(:))'));
        end   
        Params.coilNum = length(Params.coilName);
    end
    
    if verbose        
        disp(['searching DICOM ', num2str(ii), ' ...']);
    end
end

TE = single(zeros([NumEcho 1]));

% Number of Slices
if isfield(info, 'SpacingBetweenSlices')
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/info.SpacingBetweenSlices) + 1;
else
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/Params.voxSize(3)) + 1;
end

%% read in imaging data
if Manufacturer ~= 3   % for SIEMENS or PHILIPS or TOSHIBA, can read in Mag & Phase directly
    GREPhase = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREMag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
else
    % For GE data, read in Real/Imaginary then convert to Mag/Phase
    GREReal = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREImag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    
    GREPhase = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREMag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
end

for i = 1:length(filelist)
    info = dicominfo([DICOMdir '/' filelist(i).name]);

    if isfield(info, 'SpacingBetweenSlices') && Manufacturer ~= 3   % slice
        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/info.SpacingBetweenSlices) +1);      % slice number
    else
        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/Params.voxSize(3)) +1);      % slice number
    end

    if TE(info.EchoNumber)==0
        TE(info.EchoNumber)=info.EchoTime;     % in ms     % echo time
    end
    
    if isfield(info, 'TemporalPositionIdentifier')          % dynamics
        dynamic = info.TemporalPositionIdentifier;
    else
        dynamic = 1;
    end

    if isfield(info, 'Private_0051_100f')                   % SIEMENS coil
        coil = find(strcmpi(Params.coilName, cellstr(char((info.Private_0051_100f(:))')))>0, 1, 'first');        % coil name
    else
        coil = 1;
    end        

    % read in image data
    if Manufacturer < 3         % for SIEMENS and PHILIPS data
        if strcmpi(info.ImageType(18), 'P')             % phase
            switch Manufacturer
                case 1
                    GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread([DICOMdir '/' filelist(i).name])')*info.RescaleSlope+info.RescaleIntercept)/single(info.LargestImagePixelValue)*pi;    %phase
                case 2
                    if isfield(info, 'RealWorldValueMappingSequence')
                        GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = 1e-3*(single(dicomread([DICOMdir '/' filelist(i).name])')*info.RealWorldValueMappingSequence.Item_1.RealWorldValueSlope+info.RealWorldValueMappingSequence.Item_1.RealWorldValueIntercept); %phase
                    elseif isfield(info, 'RescaleSlope') && strcmpi(info.RescaleType, 'milliradials')
                        GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = 1e-3*(single(dicomread([DICOMdir '/' filelist(i).name])')*info.RescaleSlope+info.RescaleIntercept); %phase
                    else
                        error('unknown scaling ...')
                    end
            end

        elseif strcmpi(info.ImageType(18), 'M')         % magnitude
            GREMag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread([DICOMdir '/' filelist(i).name])');     %magnitude
        end
    elseif Manufacturer == 3
        % for GE data
        if mod(info.InstanceNumber,4)==0
            GREImag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread([DICOMdir '/' filelist(i).name])');   %imaginary
        elseif mod(info.InstanceNumber,4)==3
            GREReal(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread([DICOMdir '/' filelist(i).name])');   %real
        end     
    elseif Manufacturer == 4
        % for TOSHIBA
        if info.WindowCenter == 0   % phase
            GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread([DICOMdir '/' filelist(i).name])'))/single(info.WindowWidth)*2*pi;
        else
            GREMag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread([DICOMdir '/' filelist(i).name])');     %magnitude
        end
    end
    
    % Only show waitbar every 10 iterations
    if(mod(i, 10) == 0)
        multiWaitbar(textWaitbar, (i/length(filelist)));
    end
end

% Close waitbar
multiWaitbar( 'CloseAll' );

if Manufacturer == 3    % for GE data with phase shift artifacts
    temp = ifft(fftshift(fft(GREReal + 1i.*GREImag, [], 3), 3), [], 3);    
    temp = conj(temp);  % for GE data, use the conj (-phase)
    GREMag = abs(temp);
    GREPhase = angle(temp);
    clear temp
end

if Manufacturer == 4        % TOSHIBA data is flipped?
    GREMag = flip(GREMag, 3);
    GREPhase = flip(GREPhase, 3);
end

% Post-process
Params.nEchoes  = NumEcho;
Params.voxSize  = (Params.voxSize);
Params.sizeVol  = double(Params.sizeVol);
Params.fov      = (round(Params.voxSize.*Params.sizeVol));
Params.TEs      = TE./1000;
Params.nDynamics= Params.DynamicNum;

% data in 6D array, ncol, nrow, necho, ndyanmics, ncoil
