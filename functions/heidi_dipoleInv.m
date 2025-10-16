function chi = heidi_dipoleInv(tissue_phase, brain_mask, B0, TE, voxel_size, magnitude_img)
% Inputs:
% phase - 3D tissue phase images
% mag - (optional) 3D magnitude image from the GRE sequence
% mask - 3D binary brain mask.
% field_strength - Main magnetic field strength in Tesla
% TE - Echo time in seconds
% voxel_size - 1x3 vector of voxel dimensions in mm

% Output:
% QSM = Reconstructed quantitative susceptibility map in ppm.

% Require:
% Image Processing Toolbox

% Reference:
% Schweser, F., Sommer, K., Deistung, A., & Reichenbach, J. R. (2012). 
% Quantitative susceptibility mapping for investigating subtle susceptibility 
% variations in the human brain. NeuroImage, 62(3), 2083–2100. 
% https://doi.org/10.1016/j.neuroimage.2012.05.067

% Check for magnitude image
if nargin < 6
    magnitude_img = [];
    fprintf('No magnitude image provided. Skipping magnitude-based prior.\n');
end

% Physical constants
gamma = 267.513e6; % Gyromagnetic ratio for protons (rad/s/T)

% Algorithm parameters (based on the paper's findings, may need tuning)
% Sub-domain definition parameters (Fig. 4h)
t_ill = 0.09;
t_trans = 0.23 - t_ill; % t_ill + t_trans = 0.23

% A priori information thresholds (Fig. 4b-d)
t_nabla_phi = 0.00105; % rad/T/s/mm
t_delta_phi = 0.0159;  % rad/T/s/mm^2
t_nabla_m_factor = 2.40; % Multiplier for the mean magnitude gradient

% Anisotropic diffusion parameters (from paper's methods)
diff_conductance = 1; % This is the 'GradientThreshold' (K value)
diff_iterations_phase = 5;
diff_iterations_chi = 8;

% L1 minimization parameters (for simplified Step 4)
l1_iterations = 50;
l1_step_size = 0.1; % Learning rate for gradient descent

% Get matrix dimensions
matrix_size = size(tissue_phase);
if ~isequal(matrix_size, size(brain_mask))
    error('Dimensions of tissue_phase and brain_mask must match.');
end

fprintf('Matrix Size: %d x %d x %d\n', matrix_size(1), matrix_size(2), matrix_size(3));
fprintf('Voxel Size: %.2f x %.2f x %.2f mm\n', voxel_size(1), voxel_size(2), voxel_size(3));

% Convert local phase to local field perturbation (in ppm)
% f = B_delta / B0 = phase / (gamma * TE * B0)
% The susceptibility 'chi' will also be in ppm.
local_field = (tissue_phase / (gamma * TE * B0)) * 1e6; % convert to ppm
local_field = local_field .* brain_mask;

% Create k-space grid
[kx, ky, kz] = create_kspace_grid(matrix_size, voxel_size);

% Calculate Dipole Kernel in k-space
% d(k) = 1/3 - kz^2 / k^2
k2 = kx.^2 + ky.^2 + kz.^2;
dipole_kernel = 1/3 - (kz.^2) ./ k2;
dipole_kernel(k2 == 0) = 0; % Set DC component to 0 to avoid division by zero

fprintf('Step 0: Initialization complete.\n');

%% --- Step 1: Identification of Spectral Sub-domains ---

% Based on Eq. (13) and the parameters derived from Fig. 4h
abs_dipole_norm = abs(dipole_kernel);

M_ill = abs_dipole_norm < t_ill;
M_trans = (abs_dipole_norm >= t_ill) & (abs_dipole_norm < t_ill + t_trans);
M_well = abs_dipole_norm >= t_ill + t_trans;

fprintf('Step 1: Identified k-space sub-domains.\n');
fprintf('         - Well-conditioned: %d voxels\n', sum(M_well(:)));
fprintf('         - Transitional:     %d voxels\n', sum(M_trans(:)));
fprintf('         - Ill-conditioned:  %d voxels\n', sum(M_ill(:)));

%% --- Step 2: Unregularized Reconstruction of Well-conditioned Sub-domain ---
% This step calculates an initial estimate, chi_interim, by simple
% deconvolution, ignoring the ill-conditioned cone.

% FFT of the local field
field_k = fftn(local_field);

% Create a modified dipole kernel for the initial inversion
% This avoids division by very small numbers in the ill-conditioned region
dipole_inv = 1 ./ dipole_kernel;
dipole_inv(M_ill) = 0; % Set inverse to 0 on the cone

% Initial susceptibility map in k-space
chi_interim_k = field_k .* dipole_inv;

% Transform back to image space
chi_interim = real(ifftn(chi_interim_k));

fprintf('Step 2: Performed initial unregularized reconstruction (chi_interim).\n');

%% --- Step 3: Extraction of a priori Homogeneity Information ---
% This is the "Homogeneity Enabled" part of HEIDI. It creates spatial
% weighting masks based on gradients of the phase and magnitude images.

% Conversion factor kappa for thresholds
kappa = 1 / (TE * B0);

% --- 3a. Phase Gradient Prior ---
fprintf('Calculating phase gradient prior...\n');
denoised_phase = imdiffusefilt(tissue_phase, 'ConductionMethod', 'exponential', ...
    'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);
[phase_grad_x, phase_grad_y, phase_grad_z] = gradient(denoised_phase, voxel_size(1), voxel_size(2), voxel_size(3));

W_x_vphi = (abs(phase_grad_x) * kappa) < t_nabla_phi;
W_y_vphi = (abs(phase_grad_y) * kappa) < t_nabla_phi;
W_z_vphi = (abs(phase_grad_z) * kappa) < t_nabla_phi;

% --- 3b. Phase Laplacian Prior ---
fprintf('Calculating phase Laplacian prior...\n');
phase_lap = del2(denoised_phase, voxel_size(1), voxel_size(2), voxel_size(3));
denoised_lap = imdiffusefilt(abs(phase_lap), 'ConductionMethod', 'exponential', ...
    'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);

W_delta_phi = (denoised_lap * kappa) < t_delta_phi;

% --- 3c. Magnitude Gradient Prior (Optional) ---
if ~isempty(magnitude_img)
    fprintf('Calculating magnitude gradient prior...\n');
    denoised_mag = imdiffusefilt(magnitude_img, 'ConductionMethod', 'exponential', ...
        'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);
    [mag_grad_x, mag_grad_y, mag_grad_z] = gradient(denoised_mag, voxel_size(1), voxel_size(2), voxel_size(3));

    % Denoise the gradients themselves as per the paper
    mag_grad_x = imdiffusefilt(abs(mag_grad_x), 'ConductionMethod', 'exponential', 'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);
    mag_grad_y = imdiffusefilt(abs(mag_grad_y), 'ConductionMethod', 'exponential', 'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);
    mag_grad_z = imdiffusefilt(abs(mag_grad_z), 'ConductionMethod', 'exponential', 'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_phase);

    % Calculate threshold based on mean gradient
    mean_mag_grad_norm = mean(sqrt(mag_grad_x(brain_mask).^2 + mag_grad_y(brain_mask).^2 + mag_grad_z(brain_mask).^2));
    t_nabla_m = t_nabla_m_factor * mean_mag_grad_norm;

    W_x_vm = abs(mag_grad_x) < t_nabla_m;
    W_y_vm = abs(mag_grad_y) < t_nabla_m;
    W_z_vm = abs(mag_grad_z) < t_nabla_m;
else
    % If no magnitude image, the mask is all ones (no contribution)
    W_x_vm = ones(matrix_size);
    W_y_vm = ones(matrix_size);
    W_z_vm = ones(matrix_size);
end

% --- 3d. Combine Priors ---
% W_j = W_j^nabla_phi * W^delta_phi * W_j^nabla_m (Eq. 10)
W_x = W_x_vphi .* W_delta_phi .* W_x_vm;
W_y = W_y_vphi .* W_delta_phi .* W_y_vm;
W_z = W_z_vphi .* W_delta_phi .* W_z_vm;

% As per the paper, add a small weight to non-homogeneous regions
% to provide weak total-variation regularization everywhere.
W_x(~W_x) = 0.1;
W_y(~W_y) = 0.1;
W_z(~W_z) = 0.1;

fprintf('Step 3: Extracted a priori homogeneity masks (W).\n');

%% --- Step 4: Reconstruction of the Ill-conditioned Sub-domain ---
% This step solves a constrained L1-minimization problem.
% min ||chi||_WG subject to P_ill * chi = P_ill * chi_interim
% The paper uses NESTA. Here, we use a simplified iterative gradient
% descent as a demonstration. For best results, a dedicated L1 solver
% (e.g., from a toolbox like TVAL3 or Split Bregman) is recommended.

fprintf('Step 4: Reconstructing ill-conditioned domain (L1 minimization)... \n');

chi_ill_recon = chi_interim; % Start with the interim solution

% Iterative gradient descent
for iter = 1:l1_iterations
    if mod(iter, 10) == 0
        fprintf('   Iteration %d/%d\n', iter, l1_iterations);
    end

    % Calculate spatial gradients of the current chi estimate
    [chi_grad_x, chi_grad_y, chi_grad_z] = gradient(chi_ill_recon, voxel_size(1), voxel_size(2), voxel_size(3));

    % Calculate the update term (subgradient of the L1 norm)
    update_term = W_x .* sign(chi_grad_x) + W_y .* sign(chi_grad_y) + W_z .* sign(chi_grad_z);

    % We need the divergence of the update term for the gradient descent step
    % div( [wx*sign(gx), wy*sign(gy), wz*sign(gz)] )
    [div_x, ~] = gradient(update_term, voxel_size(1), voxel_size(2), voxel_size(3));

    % Gradient descent step
    chi_ill_recon = chi_ill_recon - l1_step_size * div_x;
    chi_ill_recon = chi_ill_recon .* brain_mask; % Apply mask

    % Enforce k-space constraint: keep well and trans domains fixed
    chi_k_updated = fftn(chi_ill_recon);
    chi_k_updated(M_well | M_trans) = chi_interim_k(M_well | M_trans);
    chi_ill_recon = real(ifftn(chi_k_updated));
end

fprintf('Step 4: Ill-conditioned domain reconstruction complete.\n');

%% --- Step 5: Reconstruction of the Transitional Sub-domain ---
% This step denoises the result from Step 4 to clean up the transitional
% sub-domain, which is susceptible to noise amplification.

fprintf('Step 5: Denoising transitional domain...\n');

% The paper denoises the entire map from Step 4.
% Apply a mask before denoising to focus on the brain.
chi_to_denoise = chi_ill_recon;
mean_brain_val = mean(chi_to_denoise(brain_mask));
chi_to_denoise(~brain_mask) = mean_brain_val;

chi_trans_denoised = imdiffusefilt(chi_to_denoise, 'ConductionMethod', 'exponential', ...
    'GradientThreshold', diff_conductance, 'NumberOfIterations', diff_iterations_chi);

fprintf('Step 5: Denoising complete.\n');

%% --- Final Combination ---
% Assemble the final susceptibility map from the different k-space domains.

% Get k-space of the final denoised map
chi_denoised_k = fftn(chi_trans_denoised);

% Start with the k-space from the initial reconstruction
chi_final_k = chi_interim_k;

% Replace the ill-conditioned and transitional parts with the refined,
% denoised versions from steps 4 and 5.
chi_final_k(M_ill) = chi_denoised_k(M_ill);
chi_final_k(M_trans) = chi_denoised_k(M_trans);

% Inverse FFT to get the final susceptibility map
chi = real(ifftn(chi_final_k));
chi = chi .* brain_mask; % Final masking

fprintf('--- HEIDI QSM Finished ---\n');



end