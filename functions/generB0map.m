function [freq_map,num_echoes] = generB0map(unwrapped_phase_4d, mask_brain, TE_array)

[x_dim, y_dim, z_dim, num_echoes] = size(unwrapped_phase_4d);
freq_map = zeros(x_dim, y_dim, z_dim);

% Linear fit phase vs. echo time voxel-wise
fprintf('Starting multi-echo linear fit...\n');

for ix = 1:x_dim
    for iy = 1:y_dim
        for iz = 1:z_dim
            if mask_brain(ix, iy, iz) > 0
                phase_vector = squeeze(unwrapped_phase_4d(ix, iy, iz, :));
                P = polyfit(TE_array, phase_vector', 1);
                freq_map(ix, iy, iz) = P(1);
            end
        end
    end
end

fprintf('Multi-echo linear fit completed\n');

end