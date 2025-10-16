function D = dipole_kernel(B0_dir, matrix_size, voxel_size)
% Dipole kernel function

    [kx, ky, kz] = ndgrid(...
        ifftshift(-floor(matrix_size(1)/2):ceil(matrix_size(1)/2)-1)/ (matrix_size(1)*voxel_size(1)),...
        ifftshift(-floor(matrix_size(2)/2):ceil(matrix_size(2)/2)-1)/ (matrix_size(2)*voxel_size(2)),...
        ifftshift(-floor(matrix_size(3)/2):ceil(matrix_size(3)/2)-1)/ (matrix_size(3)*voxel_size(3)));

    k2 = kx.^2 + ky.^2 + kz.^2;
    k2(k2==0) = eps; % Avoid division by zero

    % Normalize B0 direction
    B0_dir = B0_dir / norm(B0_dir);

    k_dot_B0 = kx * B0_dir(1) + ky * B0_dir(2) + kz * B0_dir(3);
    
    D = (1/3) - (k_dot_B0.^2) ./ k2;
end
