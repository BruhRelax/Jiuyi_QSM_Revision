function [tissuePhase, qsm] = singleEcho (unwraped_phase, brain_mask, magnitude, info, TE)
%% Run qsm processing single Echo script
voxel_size = info.PixelDimensions;
TransfMatrix = info.Transform.T;
brain_mask = logical(brain_mask);
magnitude = single(magnitude);
matrix_size = size(unwraped_phase);
padsize = [12,12,12];
case_num = input("Enter processed Case Number:", "s");
% padsize: Identify largest matrix dimension, then Choose the next power of two, which ensures optimal FFT performance and enough padding.

bfrMethod = input ('Enter Background Field Removal method (RESHARP/VSHARP): ', 's');

if strcmp(bfrMethod, 'RESHARP')
    disp('Processing Background Field Removal (RESHARP)...')
    tic;
    [tissuePhase, mask_ero] = resharp(unwraped_phase, brain_mask, voxel_size, 4, 1e-6, 100);
    elapsedTime = toc;
    fprintf('RESHARP completed in %.2f seconds.\n', elapsedTime);

elseif strcmp(bfrMethod, 'VSHARP')
    disp('Processing Background Field Removal (VSHARP)...')
    [tissuePhase, mask_ero] = V_SHARP(unwraped_phase, brain_mask, 'voxelsize',voxel_size, 'smvsize', 12);
end

%% Parameters for MEDI
RDF = single(tissuePhase);
iFreq = single(unwraped_phase);
N_std = single(ones(size(RDF)));
Mask = brain_mask; 
iMag = magnitude;
row_vector = TransfMatrix(1:3,1);
col_vector = TransfMatrix(1:3,2);
slice_normal = cross(row_vector, col_vector);
B0_dir = round(slice_normal / norm(slice_normal));
CF = 128;

save('RDF.mat','iFreq','RDF','Mask','iMag','N_std','matrix_size','voxel_size','CF','B0_dir');

%% Dipole inversion
diMethod = input('Enter dipole inversion method (MEDI/iLSQR): ', 's');

if strcmp(diMethod, 'MEDI')
    % Dipole Inversion (MEDI)
    disp('Processing Dipole Inversion (MEDI)...')
    run('toolboxes/MEDI_toolbox/MEDI_set_path.m');
    qsm = MEDI_L1('lambda',1000,'smv',12);

elseif strcmp(diMethod, 'iLSQR')
    % Dipole inversion (iLSQR)
    disp('Processing Dipole Inversion (iLSQR)...')
    qsm = QSM_iLSQR(tissuePhase,brain_mask,'TE', TE,'B0', 3,'H',B0_dir,'padsize',padsize,'voxelsize',voxel_size); 

elseif strcmp(diMethod, 'HEIDI')
    disp('Processing Dipole Inversion (HEIDI)...')
    qsm = heidi_dipoleInv(tissuePhase, brain_mask, 3, TE/1000, voxel_size, magnitude);
else
    disp('Invalid Method')

end

outputDir = fullfile('output_data', case_num);
mkdir(outputDir);
tissuePhaseFilename = sprintf('tissue_phase_%s.nii.gz', bfrMethod);
qsmFilename = sprintf('QSM_%s_%s.nii.gz', diMethod, bfrMethod);

% --- Construct the final, full file paths ---
tissuePhasePath = fullfile(outputDir, tissuePhaseFilename);
qsmPath = fullfile(outputDir, qsmFilename);

niftiwrite(single(tissuePhase),tissuePhasePath, info);
niftiwrite(single(qsm),qsmPath, info);

end