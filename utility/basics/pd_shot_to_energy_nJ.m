% File: pd_shot_to_energy_nJ.m
% Weiye Li - Razansky Lab - University and ETH Zurich
% E-mail: weiye.li@uzh.ch
% Date: March 25, 2023

function ppenJShots = pd_shot_to_energy_nJ(pdShots, currWavelength)
%converts 1D photodiode signal to per-pulse energy in nJ, based on calibration data against power meter

% INPUT: 
% pdShots: photodiode signals for every laser shot, can be a 2D or 3D array, first dimension is time (t)
% pdRiseTime: rise time of the photodiode specified by Thorlabs, [ns]
% dt: sampling time interval, [s]

% OUTPUT: per-pulse energy in nJ for every laser shot

  % reshape to [nt, nShots] if input is 3D
  pdShots = single(pdShots);
  nDim = numel(size(pdShots));
  if nDim == 3
      [ntPd, nx, ny] = size(pdShots);
      pdShots = reshape(pdShots, [ntPd, nx*ny]);
  end

  [~, maxIndices] = max(pdShots, [], 1);
  maxIndices = unique(maxIndices); % also corrected for laser jitter, so this should only be a single number
  if length(maxIndices) > 1
    error('Did you do laser jitter correction?');
  end

  ppeauShots = sum(pdShots, 1); % [1, nShots]
  % NOTE: taking max. does NOT correlate with power meter readings well!

  if currWavelength == 532
    pFit_532 = [0.0263, 40.3815];
    ppenJShots = polyval(pFit_532, ppeauShots);
  elseif currWavelength == 558
    pFit_558 = [0.0269, -20.5326];
    ppenJShots = polyval(pFit_558, ppeauShots);
  else
    error('Invalid wavelength!');
  end

  if nDim == 3
    ppenJShots = reshape(ppenJShots, [nx, ny]); % ppe in arbitrary unit
  end

end