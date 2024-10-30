% File: Load_Raw_Data.m @ dualwaveProcessor
% Date: July 28th, 2022
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Load_Raw_Data(dwp)

  fprintf('[dualwaveProcessor] Loading raw data...');

  D = matfile(fullfile(dwp.data_dir, dwp.data_name));
  dwp.scanSett = D.sett;

  % apply FOV cropping to get rid of edge area (e.g. first and last 0.5mm in both x and y)
  if dwp.flagCropFOVBoundary
    warning('Enabling FOV cropping...');
    scanSettTmp = D.sett;
    scanStepSize_mm = scanSettTmp.dr(1);
    assert(scanStepSize_mm == scanSettTmp.dr(2));
    nCropX = dwp.cropXBoundary_mm / scanStepSize_mm;
    nCropY = dwp.cropYBoundary_mm / scanStepSize_mm;
    nx = scanSettTmp.nALines;
    ny = scanSettTmp.nBScans;
    if nCropY == 0
      % only crop along fast stage because of PRF fluctuation
      dwp.RawDataPd = D.RawDataPd(:, :, :, (nCropX+1):(nx-nCropX), :);
      dwp.RawDataUs = D.RawDataUs(:, :, :, (nCropX+1):(nx-nCropX), :);
    else
      % can also crop along slow stage because of low SNR at FOV edges
      dwp.RawDataPd = D.RawDataPd(:, :, :, (nCropX+1):(nx-nCropX), nCropY:(ny-nCropY));
      dwp.RawDataUs = D.RawDataUs(:, :, :, (nCropX+1):(nx-nCropX), nCropY:(ny-nCropY));
    end
    clear scanSettTmp;
  else
    dwp.RawDataPd = D.RawDataPd;
    dwp.RawDataUs = D.RawDataUs;
  end

  clear D;

  fprintf('done!\n');

end