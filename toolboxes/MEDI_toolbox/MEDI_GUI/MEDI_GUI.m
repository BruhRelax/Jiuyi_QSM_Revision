function varargout = MEDI_GUI(varargin)
% MEDI_GUI MATLAB code for MEDI_GUI.fig
%      MEDI_GUI, by itself, creates a new MEDI_GUI or raises the existing
%      singleton*.
%
%      H = MEDI_GUI returns the handle to a new MEDI_GUI or the handle to
%      the existing singleton*.
%
%      MEDI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEDI_GUI.M with the given input arguments.
%
%      MEDI_GUI('Property','Value',...) creates a new MEDI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEDI_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEDI_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEDI_GUI

% Last Modified by GUIDE v2.5 22-Nov-2020 16:47:37
% Last Modified by Ramin 5/2/20

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEDI_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MEDI_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MEDI_GUI is made visible.
function MEDI_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEDI_GUI (see VARARGIN)

% Choose default command line output for MEDI_GUI
handles.output = hObject;
handles.unwrapSelection = '';
handles.varin = {};
handles.loadedDataFolder = '';
handles.filedir = '';
handles.resultFolder = '';
handles.SMVEnable = 0;
handles.currPath = '';
handles.SMVEnable = 1;
handles.MERIT = 1;
handles.lambda = 1000;
handles.percentage = 0.9;
handles.SMVradius = 5;
%handles.defaults.selectedobject.uipanel13=handles.advRB; % runmode
handles.defaults.selectedobject.uipanel6=handles.BET; % BFR
handles.defaults.selectedobject.uipanel9=handles.RGRadioButton;
% handles.defaults.selectedobject.uipanel9=handles.LaplacianButton;
handles.defaults.selectedobject.BFRSel=handles.PDFButton;
% handles.defaults.selectedobject.BFRSel=handles.LBVButton;
handles.defaults.value.SMVbox=handles.SMVEnable;
handles.defaults.value.MERITbox=handles.MERIT;
handles.defaults.string.lambdaEditText=num2str(handles.lambda);
handles.defaults.string.SMVEdit=num2str(handles.SMVradius);
handles.defaults.string.edgeEdit=num2str(handles.percentage);
handles.defaults.string.MaskLocationEditable='';

handles.runmode='advRB'; 
set(handles.uipanel13,'selectedobject',handles.(handles.runmode));
set_defaults(handles);
enable_elements(handles, handles.runmode);
% 
% % default is BET
% set(handles.uipanel6,'selectedobject',handles.BET);
% % default is SMV
% 
% set(handles.SMVbox,'value',handles.SMVEnable);
% % default is MERIT
% 
% set(handles.MERITbox,'value',handles.MERIT);
% Update handles structure
% assignin('base','handles_copy', handles); 
guidata(hObject, handles);

% UIWAIT makes MEDI_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEDI_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function set_defaults(handles)
for kind=fieldnames(handles.defaults)'
    for f=fieldnames(handles.defaults.(kind{:}))'
        set(handles.(f{:}),kind{:},handles.defaults.(kind{:}).(f{:}))
    end
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseDataButton.
function BrowseDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExPath= uigetdir();

% ExPath = fullfile(FilePath, FileName);
set(handles.DataLocationEditable,'string',ExPath);


function DataLocationEditable_Callback(hObject, eventdata, handles)
% hObject    handle to DataLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataLocationEditable as text
%        str2double(get(hObject,'String')) returns contents of DataLocationEditable as a double


% --- Executes during object creation, after setting all properties.
function DataLocationEditable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over BrowseDataButton.
function BrowseDataButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to BrowseDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = Load(handles)
handles.files = [];
FilePath = get(handles.DataLocationEditable,'string');
set(handles.statusText,'String','Busy Loading Data...  Please Wait');
drawnow
cmdstr=['[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,files_local]' ...
    ' = Read_DICOM(''' FilePath ''');'];
evalin('base', cmdstr);
handles.files = evalin('base', 'files_local');
set(handles.statusText,'String','Done Loading Data');
handles.loadedDataFolder = FilePath;
% assignin('base','iField',iField);
% assignin('base','voxel_size',voxel_size);
% assignin('base','matrix_size',matrix_size);
% assignin('base','CF',CF);
% assignin('base','delta_TE',delta_TE);
% assignin('base','TE',TE);
% assignin('base','B0_dir',B0_dir);
% assignin('base','files_local',files_local);

% --- Executes on button press in Loadbutton.
function Loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, Load(handles));

function handles = FitPPM(handles)
% iField = evalin('base','iField');
% TE = evalin('base','TE');
set(handles.statusText,'String','Busy Fitting...');
drawnow
evalin('base','[iFreq_raw, N_std] = Fit_ppm_complex(iField);');
% evalin('base','[iFreq_raw, N_std] = Fit_ppm_complex_TE(iField,TE);');
set(handles.statusText,'String','Fitting Done');
% assignin('base','iFreq_raw',iFreq_raw);
% assignin('base','N_std',N_std);

% --- Executes on button press in FitButton.
function FitButton_Callback(hObject, eventdata, handles)
% hObject    handle to FitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FitPPM(handles);

% % --- Executes on button press in pushbutton5.
% function pushbutton5_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton5 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles.Mask = genMask(handles.iField, handles.voxel_size);


% --- Executes on button press in BrowseMaskButton.
function BrowseMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileNameMask,FilePathMask]= uigetfile();
MaskFile = fullfile(FilePathMask, FileNameMask);
set(handles.MaskLocationEditable,'string',MaskFile);


function MaskLocationEditable_Callback(hObject, eventdata, handles)
% hObject    handle to MaskLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskLocationEditable as text
%        str2double(get(hObject,'String')) returns contents of MaskLocationEditable as a double


% --- Executes during object creation, after setting all properties.
function MaskLocationEditable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = Mask(handles)
type = get(get(handles.uipanel6,'SelectedObject'),'Tag');
fprintf('%s Mask\n', type);
switch type
    case 'Auto'
%         iField = evalin('base','iField');
%         voxel_size = evalin('base','voxel_size');
%         mask = type;
        evalin('base', 'Mask = genMask(iField, voxel_size);');
        assignin('base','mask',type);
%         assignin('base','Mask',Mask);
        set(handles.statusText,'String','Mask Automatically Generated');
        drawnow
    case 'BET'
%         iField = evalin('base','iField');
%         voxel_size = evalin('base','voxel_size');
%         matrix_size = evalin('base','matrix_size');
%         mask = type;
%         TE = evalin('base','TE');
%         try
%             iMag = evalin('base','iMag');
%         catch
%             iMag = sqrt(sum(abs(iField).^2,4));
%             assignin('base','iMag',iMag);
%         end
        evalin('base', 'iMag = sqrt(sum(abs(iField).^2,4));');
        evalin('base', 'Mask = BET(iMag, matrix_size, voxel_size);');
        evalin('base', 'R2s = arlo(TE, abs(iField));');
        evalin('base', 'Mask_CSF = extract_CSF(R2s, Mask, voxel_size);');
        assignin('base','mask',type);
        set(handles.statusText,'String','Mask Automatically Generated');
        drawnow
    case 'UserSelect'
        mask = type;
        assignin('base','mask',mask);
        Mask = loadMask(handles);
        if ~isempty(Mask)
            assignin('base','Mask',Mask);
        end
    otherwise
        disp(':(');
end

function handleMaskBrowse(handles, type)
switch type
    case 'Auto'
        set(handles.MaskLocationEditable,'Enable','off');
        set(handles.BrowseMaskButton,'Enable','off');
        set(handles.LoadMaskButton,'Enable','on');
        drawnow
    case 'BET'
        set(handles.MaskLocationEditable,'Enable','off');
        set(handles.BrowseMaskButton,'Enable','off');
        set(handles.LoadMaskButton,'Enable','on');
        drawnow
    case 'UserSelect'
        set(handles.MaskLocationEditable,'Enable','on');
        set(handles.BrowseMaskButton,'Enable','on');
        set(handles.LoadMaskButton,'Enable','on');
    otherwise
        disp(':(');
end

% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
handleMaskBrowse(handles, get(eventdata.NewValue,'Tag'));

function handles = Unwrap(handles)
type = get(get(handles.uipanel9,'SelectedObject'),'Tag');
fprintf('%s Unwrap\n', strrep(type, 'RadioButton', ''));
switch type
    case 'LaplacianRadioButton'
        set(handles.statusText,'String','Busy Performing Laplacian Unwrapping');
        drawnow
        evalin('base', 'iFreq = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size);');
        set(handles.statusText,'String','Done Performing Laplacian Unwrapping');
    case 'RGRadioButton'
        set(handles.statusText,'String','Busy Performing Region Growth Unwrapping');
        drawnow
        evalin('base', 'iFreq = unwrapPhase(iMag,iFreq_raw,matrix_size);');
        set(handles.statusText,'String','Done Performing Region Growth Unwrapping');
    otherwise 
        disp('fail');
end

% --- Executes on button press in UnwrapButton.
function UnwrapButton_Callback(hObject, eventdata, handles)
% Laplacian
% hObject    handle to UnwrapButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    Unwrap(handles);
catch
    FitButton_Callback(hObject, eventdata, handles);
    UnwrapButton_Callback(hObject, eventdata, handles);
end

function handles=RemoveBackground(handles)
type = get(get(handles.BFRSel,'SelectedObject'),'Tag');
fprintf('Remove background using %s\n', strrep(type, 'Button', ''));
switch type
    case 'PDFButton'
        set(handles.statusText,'String','Busy Performing BFR...');
        drawnow
        evalin('base', 'RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);');
        set(handles.statusText,'String','Done Performing BFR');
%         catch ME2
%             disp(['ME2' ME2.identifier]);
%             UnwrapButton_Callback(hObject, eventdata, handles);
%             RDFButton_Callback(hObject, eventdata, handles)
%         end
    case 'LBVButton'
%         try
        set(handles.statusText,'String','Busy Performing BFR...');
        drawnow
        evalin('base', 'RDF = LBV(iFreq,Mask,matrix_size,voxel_size);');
        set(handles.statusText,'String','Done Performing BFR');
%         catch ME3
%             disp(['ME3' ME3.identifier]);
%             UnwrapButton_Callback(hObject, eventdata, handles);
%             RDFButton_Callback(hObject, eventdata, handles);
%         end
    otherwise
        disp('Error');
end

% --- Executes on button press in PDFButton.
function RDFButton_Callback(hObject, eventdata, handles)
% hObject    handle to PDFButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RemoveBackground(handles);

function handles = MEDI(handles, args)
if nargin<2
    args={};
end
handles.varin={};
tbVal = @(x) str2num(get(handles.(x),'String'));
tbValid = @(x) all(ismember(get(handles.(x),'String'), '0123456789+-.eEdD'));
if tbValid('SMVEdit')
    if(handles.SMVEnable)
        handles.varin = {'smv', tbVal('SMVEdit'), handles.varin{:}};
        disp('SMV Entered');
    else
        handles.varin = {'nosmv', handles.varin{:}};
        disp('SMV Disabled');
    end
    %     disp(handles.varin);
end
if(handles.MERIT)
    handles.varin = {'merit'};
    disp('MERIT Enabled');
else
    handles.varin = {'nomerit'};
    disp('MERIT Disabled');
end
if tbValid('lambdaEditText')
    handles.varin = {'lambda', tbVal('lambdaEditText'), handles.varin{:}};
end
if tbValid('edgeEdit')
    handles.varin = {'percentage', tbVal('edgeEdit'), handles.varin{:}};
end
try
    set(handles.statusText,'String','Busy Performing MEDI...');
    drawnow
    QSM = MEDI_L1('filename', handles.filedir, ...
        'resultsdir', fullfile(handles.matFolder, 'results'),...
        handles.varin{:}, args{:});
    set(handles.statusText,'String','Done Performing MEDI');
    assignin('base','QSM',QSM);
catch ME
     disp(ME.identifier);
     handles = Save(handles);
     handles = MEDI(handles, args);
end

% --- Executes on button press in MEDIButton.
function MEDIButton_Callback(hObject, eventdata, handles)
% hObject    handle to MEDIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handles.varin = {handles.varin{:}, 'filename', handles.filedir,...
%     'resultsdir', fullfile(handles.matFolder, 'results')};
handles = MEDI(handles);
guidata(hObject, handles);
    


% --- Executes on button press in Visu3DButton.
function Visu3DButton_Callback(hObject, eventdata, handles)
% hObject    handle to Visu3DButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QSM = evalin('base','QSM');
voxel_size = evalin('base','voxel_size');
Mask = evalin('base','Mask');
Visu3D( QSM.*Mask, 'dimension',voxel_size);


% --- Executes when selected object is changed in uipanel9.
function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel9 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'LaplacianRadioButton'
        handles.unwrapSelection = get(eventdata.NewValue,'Tag');
    case 'RGRadioButton'
        handles.unwrapSelection = get(eventdata.NewValue,'Tag');
    otherwise
        disp(':(');
end


% --- Executes during object creation, after setting all properties.
function uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function lambdaEditText_Callback(hObject, eventdata, handles)
% hObject    handle to lambdaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambdaEditText as text
%        str2double(get(hObject,'String')) returns contents of lambdaEditText as a double


% --- Executes during object creation, after setting all properties.
function lambdaEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambdaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edgeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to edgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edgeEdit as text
%        str2double(get(hObject,'String')) returns contents of edgeEdit as a double


% --- Executes during object creation, after setting all properties.
function edgeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SMVEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SMVEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SMVEdit as text
%        str2double(get(hObject,'String')) returns contents of SMVEdit as a double


% --- Executes during object creation, after setting all properties.
function SMVEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SMVEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function enable_elements(handles, type)
hlist={'FitButton','UnwrapButton','RDFButton','MEDIButton',...
    'LaplacianRadioButton','RGRadioButton','PDFButton',...
    'LBVButton','SMVEdit','lambdaEditText','edgeEdit',...
    'saveRDFButton','SMVbox', 'MERITbox', 'LoadMaskButton',...
    'BrowseMaskButton', 'StripButton','MaskLocationEditable',...
    'Auto', 'BET', 'UserSelect'};
switch type
    case 'defRB'
        cellfun(@(x)set(handles.(x),'Enable','off'), hlist);
        set(handles.buttonJust,'Enable','on');
        set_defaults(handles);
    case 'advRB'
        cellfun(@(x)set(handles.(x),'Enable','on'), hlist);
        handleMaskBrowse(handles, get(get(handles.uipanel6,'SelectedObject'),'Tag'));
        set(handles.buttonJust,'Enable','off');
    otherwise
        disp(':(');
end


% --- Executes when selected object is changed in uipanel13.
function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel13 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
enable_elements(handles,get(eventdata.NewValue,'Tag'))


function handles = pipeline(handles)
handles = Load(handles);
handles = Mask(handles);
handles = FitPPM(handles);
handles = Unwrap(handles);
handles = RemoveBackground(handles);
handles = Save(handles);
handles = MEDI(handles);

% --- Executes on button press in buttonJust.
function buttonJust_Callback(hObject, eventdata, handles)
% hObject    handle to buttonJust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clear -excetp handles.DataLocationEditable handles.statusText
guidata(hObject, pipeline(handles));


function Mask = loadMask(handles) 
FilePath = get(handles.MaskLocationEditable,'string');
[pathstr,name,ext] = fileparts(FilePath);

switch ext
    case '.img'
        disp('Loading Mask');
        matrix_size = evalin('base','matrix_size');
        fid = fopen(FilePath);
        Mask = fread(fid,inf,'ushort');
        fclose(fid);
        Mask = reshape(Mask,matrix_size);
        disp('Done');
    case '.mat'
        disp('Loading Mask');
        mmask = load(FilePath,'-mat');
        Mask = mmask.Mask;
        disp('Done');
    otherwise 
        Mask = [];
        disp('invalid');
end

% --- Executes on button press in LoadMaskButton.
function LoadMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mask(handles);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateLBButton.
function updateLBButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateLBButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
set(handles.listbox1,'String',vars);


% --- Executes on button press in visButton.
function visButton_Callback(hObject, eventdata, handles)
% hObject    handle to visButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varName = get(handles.listbox1,'String');
varValue = get(handles.listbox1,'Value');
varrrr = varName(varValue);
varrr = cell2mat(varrrr);
afk = evalin('base',varrr);
vis(afk);

function handles = makeResultFolder(handles)
handles.resultFolder = [get(handles.DataLocationEditable,'string') '_result'];
if ~exist(handles.resultFolder,'dir')
    mkdir(handles.resultFolder);
end

function handles = makeMatFolder(handles)
handles.matFolder = [get(handles.DataLocationEditable,'string') '_mat'];
if ~exist(handles.matFolder,'dir')
    mkdir(handles.matFolder);
end

function handles=Save(handles)
handles = makeMatFolder(handles);
handles.filedir = fullfile(handles.matFolder, 'RDF.mat');
cmdstr=['save ' handles.filedir ' RDF iFreq iFreq_raw iMag N_std Mask matrix_size ' ...
        'voxel_size delta_TE CF B0_dir'];
if ismember('Mask_CSF',evalin('base','who'))
    cmdstr=[ cmdstr ' Mask_CSF'];
end
evalin('base', cmdstr);

% --- Executes on button press in saveRDFButton.
function saveRDFButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveRDFButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
guidata(hObject, Save(handles));


function handles=Write(handles)
handles = makeResultFolder(handles);
set(handles.statusText,'String',['     saving to ' handles.resultFolder ]);
drawnow
evalin('base', ['Write_DICOM(QSM, files_local,''' handles.resultFolder ''' );']);
set(handles.statusText,'String',['Done saving to ' handles.resultFolder ]);
drawnow


% --- Executes on button press in saveDICOMButton.
function saveDICOMButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveDICOMButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, Write(handles));


% --- Executes on button press in SMVbox.

function SMVbox_Callback(hObject, eventdata, handles)
% hObject    handle to SMVbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SMVbox
handles.SMVEnable = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function statusText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in StripButton.
function StripButton_Callback(hObject, eventdata, handles)
% hObject    handle to StripButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
evalin('base', 'Mask = stripBD(Mask,1);');
set(handles.statusText,'String','Mask Stripped');



function MultiEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MultiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MultiEdit as text
%        str2double(get(hObject,'String')) returns contents of MultiEdit as a double


% --- Executes during object creation, after setting all properties.
function MultiEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MultiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MultiBrowse.
function MultiBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to MultiBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paths= uigetdir();

if(isempty(get(handles.listbox2,'String')))
    new_name = [{paths}];
    set(handles.listbox2,'String',new_name);
    assignin('base','new_name',new_name);
else
    initial_name = cellstr(get(handles.listbox2,'String'));
    new_name = [initial_name;{paths}];
    set(handles.listbox2,'String',new_name);
    assignin('base','new_name',new_name);
end


% --- Executes on button press in MultiProcess.
function MultiProcess_Callback(hObject, eventdata, handles)
% hObject    handle to MultiProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear RDF.mat
new_name = evalin('base','new_name');
[x,y]=size(new_name);
for i = 1:x
    X = sprintf('%d of %d is being processed.',i,x);
    disp(X);
    FilePath = new_name{i,1};
    set(handles.DataLocationEditable,'string',FilePath);
    handles = pipeline(handles);
    handles = Write(handles);
    guidata(hObject, handles);
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearDataDirButton.
function clearDataDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearDataDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox2, 'String', '');
evalin('base', 'clearvars new_name');


% --- Executes during object creation, after setting all properties.
function uipanel13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in MERITbox.
function MERITbox_Callback(hObject, eventdata, handles)
% hObject    handle to MERITbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MERITbox
handles.MERIT = get(hObject,'Value');
guidata(hObject, handles);
