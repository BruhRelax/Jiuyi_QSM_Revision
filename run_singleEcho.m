clc;
clear;

%% load necessary Toolboxes
addpath(genpath('toolboxes'));
% addpath(genpath('toolboxes/MEDI_toolbox'));
addpath(genpath('toolboxes/QSM_SunTB'));
addpath(genpath('toolboxes/STISuite_V3.0'));
addpath(genpath('functions'));


unwraped_phase = niftiread("output_data/mse{caseID}/{unwrapped_phase}");
brain_mask = niftiread('input_data/mse{caseID}/brain_mask.nii.gz');
magnitude = niftiread('input_data/mse{caseID}/mag_robust.nii.gz');

% recommend to use header from unwrapped image
header = niftiinfo("output_data/mse{caseID}/{unwrapped_phase}");
TE = checkTEs('correspounding dicom files'); % TE unit in ms.
[tissuePhase, qsm] = singleEcho(unwraped_phase, brain_mask, magnitude, header, TE);
