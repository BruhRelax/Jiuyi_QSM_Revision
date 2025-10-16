function [kx, ky, kz] = create_kspace_grid(matrix_size, voxel_size)
    % Creates k-space coordinates for a given matrix size and voxel size.
    N = matrix_size;
    res = voxel_size;
    
    [ix, iy, iz] = ndgrid( (0:N(1)-1) - floor(N(1)/2), ...
                           (0:N(2)-1) - floor(N(2)/2), ...
                           (0:N(3)-1) - floor(N(3)/2) );

    kx = (2*pi/res(1)) * (ix/N(1));
    ky = (2*pi/res(2)) * (iy/N(2));
    kz = (2*pi/res(3)) * (iz/N(3));
end