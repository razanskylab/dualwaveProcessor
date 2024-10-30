% File: Load_PSF.m @ dualwaveProcessor
% Date: August 23, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Load_PSF(dwp)

  fprintf('[dualwaveProcessor] Loading measured point spread function...');

  psfFile = matfile(dwp.psfFile_dir);
  
  dwp.psf532 = psfFile.sphereModel532;
  dwp.psfDye = psfFile.sphereModel558;

  dwp.psfVoxel_mm = psfFile.zVoxel_mm;
  dwp.psfVec_mm_532 = psfFile.zVec_mm_532;
  dwp.psfVec_mm_dye = psfFile.zVec_mm_558;

  fprintf('done!\n');

end