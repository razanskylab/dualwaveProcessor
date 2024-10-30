function phi_vec = estim_dermis_attenuation(dermisThickness_mm, miu_s)
%Estimates fluence attenuation profile through the dermis

% INPUT: - dermisThickness_mm: assumed dermis thickness in [mm]
%        - miu_s: assumed reduced scattering coefficient in [cm-1]
% OUTPUT: phi_vec: the fluence attenuation profile in the dermis

miu_eff = miu_s; % effective attenuation in the dermis is dominated by scattering

zVoxelSize_mm = 1e-6 * 1e3; % use 1 micrometer voxel size to calculate fluence profile for good sampling
nVoxel = dermisThickness_mm / zVoxelSize_mm;
miu_eff_vec = ones(size(0:nVoxel-1)) .* miu_eff;
phi_vec = exp(-(cumsum(miu_eff_vec).*(zVoxelSize_mm*1e-1)));

end