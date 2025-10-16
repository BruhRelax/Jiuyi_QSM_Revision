function cal_b0dir_multiEcho
% getB0dirFromDICOM calculates the B0 direction vector from a DICOM file
%
% Input:
%   dicom_file: full path to a single DICOM file
%
% Output:
%   B0_direction: unit vector indicating B0 field direction

% Extract DICOM header information
info = dicominfo(dicom_file);

% Image Orientation Patient
IOP = info.ImageOrientationPatient;

% Extract row and column direction vectors
row_vector = IOP(1:3);
col_vector = IOP(4:6);

% Compute slice normal (cross product of row and column vectors)
slice_normal = cross(row_vector, col_vector);

% Normalize the slice normal to obtain B0 direction
B0_direction = slice_normal / norm(slice_normal);

end

end