function phi_vec = estim_epidermis_attenuation(epidermisThickness_mm, miu_a, miu_s)
%Estimates fluence attenuation profile through the epidermis

% INPUT: - epidermisThickness_mm: assumed epidermis thickness in [mm]
%        - miu_a: assumed absorption coefficient in [cm-1]
%        - miu_s: assumed reduced scattering coefficient in [cm-1]
% OUTPUT: phi_vec: the fluence attenuation profile in the epidermis

miu_eff = sqrt(3 * miu_a * (miu_s+miu_a)); % effective attenuation coeff.

zVoxelSize_mm = 1e-6 * 1e3; % use 1 micrometer voxel size to calculate fluence profile for good sampling
nVoxel = epidermisThickness_mm / zVoxelSize_mm;
miu_eff_vec = ones(size(0:nVoxel-1)) .* miu_eff;
phi_vec = exp(-(cumsum(miu_eff_vec).*(zVoxelSize_mm*1e-1)));

end