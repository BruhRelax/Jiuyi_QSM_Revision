function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir,files]=Read_DICOM_Siemens_3D(DicomFolder, filelist)

files=struct;

% nfiles=length(filelist);
% idx=zeros(nfiles,1);
% for i=1:nfiles
%     idx(i) = ~filelist(i).isdir;
% end
% filelist=filelist(idx>0);
nfiles=length(filelist);

filename=fullfile(DicomFolder,filelist(1).name);
info=dicominfo(filename);

item{1}.filename=filename;
files.InPlanePhaseEncodingDirection=info.SharedFunctionalGroupsSequence.Item_1.MRFOVGeometrySequence.Item_1.InPlanePhaseEncodingDirection;
files.Manufacturer = info.Manufacturer;
files.InstitutionName = info.InstitutionName;

% if ~contains(info.Manufacturer,'philips','IgnoreCase',true)
if ~contains(lower(info.Manufacturer), 'siemens')
    error('This is not a Siemens DICOM file')
end

matrix_size(1) = single(info.Width);
matrix_size(2) = single(info.Height);
voxel_size(1,1) = single(info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1));
voxel_size(2,1) = single(info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2));
voxel_size(3,1)=info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
CF = info.SharedFunctionalGroupsSequence.Item_1.MRImagingModifierSequence.Item_1.TransmitterFrequency *1e6;
NumEcho = info.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1.EchoTrainLength;

minSlice = 1e10;
maxSlice = -1e10;
f=fieldnames(info.PerFrameFunctionalGroupsSequence);
for i = 1:length(f)
    SliceLocation = info.PerFrameFunctionalGroupsSequence.(f{i}).FrameContentSequence.Item_1.InStackPositionNumber;
    ImagePositionPatient=info.PerFrameFunctionalGroupsSequence.(f{i}).PlanePositionSequence.Item_1.ImagePositionPatient;
    if SliceLocation<minSlice
       minSlice = SliceLocation;
       minLoc = ImagePositionPatient;
    end
    if SliceLocation>maxSlice
       maxSlice = SliceLocation;
       maxLoc = ImagePositionPatient;
    end
end
matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;

Affine2D = reshape(info.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient,[3 2]);
Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
B0_dir = Affine3D\[0 0 1]';
files.Affine3D=Affine3D;
files.minLoc=minLoc;
files.maxLoc=maxLoc;

itemdata = zeros([matrix_size 2*NumEcho],'single');
item=cell(2*NumEcho,1);

for j=1:nfiles
    % we already have this for the first one
    if j>1
        filename=fullfile(DicomFolder,filelist(j).name);
        item{j}.filename=filename;
        info=dicominfo(filename);
        f=fieldnames(info.PerFrameFunctionalGroupsSequence);
    end
    item{j}.EchoTime=...
        info.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime;
    item{j}.ImageType = upper(info.ImageType);
        
    fprintf('Reading pixel data from %s...\n', filename);
    data=single(dicomread(filename));
%     transpose X, Y to be consistent with Read_DICOM
    permarg=1:numel(size(data)); permarg(1:2)=[2 1];
    data=permute(data,permarg);
    
    for i = 1:length(f)
        ImagePositionPatient=info.PerFrameFunctionalGroupsSequence.(f{i}).PlanePositionSequence.Item_1.ImagePositionPatient;
        slice = int32(round(norm(ImagePositionPatient-minLoc)/voxel_size(3)) +1);
        slope = info.PerFrameFunctionalGroupsSequence.(f{i}).PixelValueTransformationSequence.Item_1.RescaleSlope;
        intercept = info.PerFrameFunctionalGroupsSequence.(f{i}).PixelValueTransformationSequence.Item_1.RescaleIntercept;
        item{j}.slice2index{slice}=f{i};
        if item{j}.ImageType(18)=='P'
            itemdata(:,:,slice,j) = (data(:,:,1,i)*slope+intercept)/abs(intercept)*pi;%phase
        elseif item{j}.ImageType(18)=='M'
            itemdata(:,:,slice,j) = data(:,:,1,i);%magnitude
        end
    end
    clear('data');
end
files.phasesign = -1;
files.zchop = 0;

EchoTimes=cellfun(@(x)x.EchoTime,item);
TE=sort(unique(EchoTimes));
if length(TE)==1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end

iField=ones([matrix_size NumEcho],'single');
for j=1:nfiles
    echo=find(TE==item{j}.EchoTime);
    if item{j}.ImageType(18)=='P'
        iField(:,:,:,echo) = iField(:,:,:,echo).*...
            exp(complex(0,-itemdata(:,:,:,j)));
    elseif item{j}.ImageType(18)=='M'
        iField(:,:,:,echo) = iField(:,:,:,echo).*...
            itemdata(:,:,:,j);
        if echo==NumEcho
            files.info3D=item{j};
        end
    end
end
clear('itemdata','item');


if 1==mod(matrix_size(3),2)
    files.slices_added=1;
    warning('Adding empty slice at bottom of volume');
    iField=padarray(iField, [0 0 1 0], 0, 'post');
    matrix_size=matrix_size+[0 0 1];
end
if ~isstruct(DicomFolder)
    disp('SIEMENS 3D READ');
end


