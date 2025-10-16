clc;
clear;

%% load necessary Toolboxes
addpath(genpath('toolboxes'));
% addpath(genpath('toolboxes/MEDI_toolbox'));
addpath(genpath('toolboxes/QSM_SunTB'));
addpath(genpath('toolboxes/STISuite_V3.0'));
addpath(genpath('functions'));


unwraped_phase = niftiread("output_data/mse8609/unwrap_PRELUDE/prelude_Holly_registered_rigid.nii.gz");
brain_mask = niftiread('input_data/8609/brain_mask.nii.gz');
magnitude = niftiread('input_data/8609/mag_robust.nii.gz');

% recommend to use header from unwrapped image
header = niftiinfo("output_data/mse8609/unwrap_PRELUDE/prelude_Holly_registered_rigid.nii.gz");
% TE = checkTEs('/Users/jiuyizhang/Desktop/UCSF/LABs/Roland_Lab/PRL_projects/simones_failed_prls/FLAIR/mse7145/55'); % TE unit in ms.
TE = 37;
[tissuePhase, qsm] = singleEcho(unwraped_phase, brain_mask, magnitude, header, TE);