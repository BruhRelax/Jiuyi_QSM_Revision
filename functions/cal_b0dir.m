function B0_direction = cal_b0dir(unwrapped_phase, mask, TE, voxelSize)
% calculate_b0_direction estimates the B0 direction from unwrapped phase data
%
% Inputs:
%   - unwrapped_phase: 3D array of unwrapped phase images (radians)
%   - mask: binary brain mask (3D logical array)
%   - TE: echo time in seconds
%   - voxelSize: voxel dimensions [x, y, z] in mm
%
% Output:
%   - B0_direction: estimated B0 direction as a unit vector [x, y, z]

% Convert phase to frequency shift
 delta_f = double(unwrapped_phase) ./ (2*pi*TE);
 mask = logical(mask);

% Optimization settings
initial_guess = [0, 0];
lb = [0, 0];
ub = [pi, 2*pi];

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');

% Optimize B0 direction
optimal_angles = fmincon(@(angles) b0_cost(angles, delta_f, mask, voxelSize), ...
                         initial_guess, [], [], [], [], lb, ub, [], options);

% Convert optimal angles to Cartesian coordinates
theta_opt = optimal_angles(1);
phi_opt = optimal_angles(2);
B0_direction = [sin(theta_opt)*cos(phi_opt), sin(theta_opt)*sin(phi_opt), cos(theta_opt)];

end

%% Cost function
function cost = b0_cost(B0_angles, delta_f, mask, voxelSize)
    theta = B0_angles(1);
    phi = B0_angles(2);
    B0_dir = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    
    D = dipole_kernel(B0_dir, size(delta_f), voxelSize);
    
    delta_f_k = fftn(delta_f);
    field_pred = real(ifftn(delta_f_k .* D));
    
    diff = delta_f(mask) - field_pred(mask);
    
    cost = sum(diff(:).^2);
end