%% Constant Parameters for QSM_Toolbox_v3.0
% Authors: Jiri van Bergen and Xu Li
% Updated by Xu Li, 2019-06-20

handles.Params.QSMdir           = fileparts(mfilename('fullpath'));
handles.Params.QSMSettingsFile  = fullfile(handles.Params.QSMdir, '/QSM_ConstantsSaved.mat');

% check FSL installation
fsldir = getenv('FSLDIR');
if isempty(fsldir)
    disp('FSLDIR set to default: /usr/local/fsl, check fsl installation...')
    setenv('FSLDIR','/usr/local/fsl');   
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    fsldir = getenv('FSLDIR');
end        

% Check if file exists
if(exist(handles.Params.QSMSettingsFile, 'file') == 2 )
    % Advice message
    disp('Loading previouly saved parameters ...')
    disp('FOR DEFAULT PARAMETERS REMOVE THE  QSM_ConstantsSaved.mat  FILE')
    
    % Load Params, only load those defined by user
    load(handles.Params.QSMSettingsFile);            
    handles.Params.saveOutput       = Params.saveOutput;
    handles.Params.B0               = Params.B0;
    
    handles.Params.FSLFolder        = Params.FSLFolder;
    handles.Params.FSLThreshold     = Params.FSLThreshold;
    
    handles.Params.EchoAvg          = Params.EchoAvg;
    handles.Params.ErodeRadius      = Params.ErodeRadius;   % mm
    handles.Params.SHARPradius      = Params.SHARPradius;   % mm 
    handles.Params.MaskThreshold    = Params.MaskThreshold;
else
    % Basic params for GUI - default or User defined
    handles.Params.saveOutput       = 1;
    handles.Params.B0               = '3';  % Tesla    
    handles.Params.FSLFolder        = fullfile(fsldir, 'bin/'); 
    handles.Params.FSLThreshold     = 0.4;  % default BET threshold
    
    handles.Params.EchoAvg          = 1;    % defualt using EchoAvg
    handles.Params.ErodeRadius      = 1;    % mm, mask erosion
    handles.Params.SHARPradius      = 8;    % mm, SHARP kernal radius
    handles.Params.MaskThreshold    = 65;   % backup masking code
end

%% Parameters not defined by Users
% Echo limits for EchoAvg
handles.TELowerLimit_1p5T       = 0/1000; % s
handles.TEUpperLimit_1p5T       = 120/1000; % s
handles.TELowerLimit_3T         = 0/1000; % s 
handles.TEUpperLimit_3T         = 60/1000; % s
handles.TELowerLimit_7T         = 0/1000; % s
handles.TEUpperLimit_7T         = 40/1000; % s
handles.TELowerLimit_11p7T      = 0/1000; % s
handles.TEUpperLimit_11p7T      = 30/1000; % s   

% Other Params
handles.Params.thresh_tsvd      = 0.05;             % good trade-off
handles.Params.gamma            = 42.57747892e6;    % Hz/T, updated
handles.Params.fileTypes        = 'method;*.par;*.PAR;*.DIC;*.dic;*.IMA;*.ima;*.DICOM;*.dicom;*.dcm;*.mat;*.1';

% Strings
handles.textReadyLoad           = 'Ready to be reconstructed';
handles.textReconstructing      = 'Reconstructing...';

% Save
guidata(hObject, handles);