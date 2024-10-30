% File: Reconstruct_Vessel_1D.m @ dualwaveProcessor
% Date: August 23, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Reconstruct_Vessel_1D(dwp)

  fprintf('[dualwaveProcessor] Performing vessel detection and reconstruction along depth...');

  usVol532 = dwp.usVol_532;
  usVol558 = dwp.usVol_dye;

  zVoxel_mm = dwp.zVoxel_mm;
  psfVoxel_mm = dwp.psfVoxel_mm;
  zVoxelSize_mm = dwp.zVoxelSizeRecon * 1e3;
  if (zVoxelSize_mm ~= mean(diff(psfVoxel_mm)))
    error('Defined reconstruction voxel size does NOT match PSF measurement!');
  end

  psf532 = dwp.psf532;
  psf558 = dwp.psfDye;

  zVec_mm_532 = dwp.zVec_mm(1,:);
  zVec_mm_558 = dwp.zVec_mm(2,:);
  psfVec_mm_532 = dwp.psfVec_mm_532;
  psfVec_mm_558 = dwp.psfVec_mm_dye;
  nt = length(zVec_mm_532);
  if (nt ~= length(psfVec_mm_532)) || (nt ~= length(psfVec_mm_558)) || (nt ~= length(zVec_mm_532))
    error('Length of time vector dooes NOT match!');
  end

  nz = length(zVoxel_mm);
  model532 = zeros(nt, nz);
  model558 = zeros(nt, nz);

  % PSF z voxel grid may be not all the reconstruction depth level we want
  % In this case, we extrapolate
  if (psfVoxel_mm(1)>zVoxel_mm(1) && psfVoxel_mm(end)<zVoxel_mm(end))
    % assign already measured PSF
    [~,psfIdx1] = min(abs(zVoxel_mm - psfVoxel_mm(1)));
    [~,psfIdxEnd] = min(abs(zVoxel_mm - psfVoxel_mm(end)));
    model532(:,psfIdx1:psfIdxEnd) = psf532;
    model558(:,psfIdx1:psfIdxEnd) = psf558;
    % extrapolate above focus
    for iAboveFocus = 1 : (psfIdx1-1)
      model532(:,iAboveFocus) = circshift(model532(:,psfIdx1), (psfIdx1-iAboveFocus)*(-1));
      model558(:,iAboveFocus) = circshift(model558(:,psfIdx1), (psfIdx1-iAboveFocus)*(-1));
    end
    % extrapolate below focus
    for iBelowFocus = (psfIdxEnd+1) : nz
      model532(:,iBelowFocus) = circshift(model532(:,psfIdxEnd), (iBelowFocus-psfIdxEnd));
      model558(:,iBelowFocus) = circshift(model558(:,psfIdxEnd), (iBelowFocus-psfIdxEnd));
    end
  end

  if (psfVoxel_mm(1)<zVoxel_mm(1)) && (psfVoxel_mm(end)>zVoxel_mm(end))
    % PSF voxel grid covers more depth than defined reconstruction grid
    % crop the PSF grid to match the desired reconstruction grid
    % TODO
    error('PSF matrix cropping is NOT yet implemented!');
  end

  [~, nAlines_oa, nBscan_oa] = size(usVol532);
  usVol532 = reshape(usVol532, [nt, nAlines_oa*nBscan_oa]);
  usVol558 = reshape(usVol558, [nt, nAlines_oa*nBscan_oa]);
  usVolRecon532 = model532' * usVol532;
  usVolRecon558 = model558' * usVol558;
  % usVolRecon532 = envelope(usVolRecon532);
  % usVolRecon558 = envelope(usVolRecon558);
  usVolRecon532 = reshape(usVolRecon532, [nz, nAlines_oa, nBscan_oa]);
  usVolRecon558 = reshape(usVolRecon558, [nz, nAlines_oa, nBscan_oa]);

  fprintf('done!\n');

end