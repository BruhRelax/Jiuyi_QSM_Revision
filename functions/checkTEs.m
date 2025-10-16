function echo_times = checkTEs(dicom_path)
% Extract unique echo times from DICOM files and additive offset from NIfTI
%
% Inputs:
% - dicom_path: folder path containing DICOM files (*.DCM)
% - nifti_path: path to the NIfTI file (*.nii.gz)
%
% Output:
% - echo_times: unique echo times (ms)

% Get list of DICOM files
dicomFiles = dir(fullfile(dicom_path, '*.dcm'));
numFiles = length(dicomFiles);
echo_times = zeros(1, numFiles);

% Extract echo times from DICOM headers
for idx = 1:numFiles
    currentFile = fullfile(dicomFiles(idx).folder, dicomFiles(idx).name);
    info = dicominfo(currentFile);
    echo_times(idx) = info.EchoTime;
end

% Get unique echo times
echo_times = unique(echo_times);

end

