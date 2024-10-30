% File: pd_shot_to_energy_au.m
% Weiye Li - Razansky Lab - University and ETH Zurich
% E-mail: weiye.li@uzh.ch
% Date: August 18, 2022

% Change log on Feb. 19, 2023: disabled pdRiseTime and dt as inputs
% function ppeauShots = pd_shot_to_energy_au(pdShots, pdRiseTime, dt)

function ppeauShots = pd_shot_to_energy_au(pdShots)
%converts 1D photodiode signal to per-pulse energy in arbitrary unit, essentially taking the sum

% INPUT: 
% pdShots: photodiode signals for every laser shot, can be a 2D or 3D array, first dimension is time (t)
% pdRiseTime: rise time of the photodiode specified by Thorlabs, [ns]
% dt: sampling time interval, [s]

% OUTPUT: per-pulse energy in arbitrary unit for every laser shot

  % reshape to [nt, nShots] if input is 3D
  pdShots = single(pdShots);
  nDim = numel(size(pdShots));
  if nDim == 3
      [ntPd, nx, ny] = size(pdShots);
      pdShots = reshape(pdShots, [ntPd, nx*ny]);
  end

  % % identify the minimum PD signal window to get rid of wavelength-dependent noise as much as possible
  % pdMeanShot = mean(pdShots, 2); % mean PD curve across all shots
  % [~, pdMaxIdx] = max(pdMeanShot, [], 1);
  % % start index of useful PD signal
  % deltaT = dt * 1e9; % [ns]
  % pdStartIdx = pdMaxIdx - ceil(pdRiseTime/deltaT);
  % % end index of useful PD signal: check at which index the PD signal relaxes to noise floor
  % pdMeanShot_noise = pdMeanShot(pdMaxIdx+50 : end); % PD signal will relax to noise floor 50 samples after the max. index
  % pdNoiseAvg = mean(pdMeanShot_noise);
  % pdNoiseStd = std(pdMeanShot_noise);
  % pdNoiseFloor = abs(pdNoiseAvg) + 2*abs(pdNoiseStd);
  % for pdEndIdx = pdMaxIdx+1 : pdMaxIdx+50
  %   if abs(pdMeanShot(pdEndIdx)) > pdNoiseFloor
  %     continue;
  %   else
  %     break;
  %   end
  % end

  [~, maxIndices] = max(pdShots, [], 1);
  maxIndices = unique(maxIndices); % also corrected for laser jitter, so this should only be a single number
  if length(maxIndices) > 1
    warning('More than one PD max. peak detected after jitter correction!');
  end
  
%   % disable negative masking
%   pdMask = (pdShots < 0);
%   pdMask = pdMask .* repmat(((1:size(pdShots,1))'>maxIndices), 1, size(pdShots,2));
%   pdMask = logical(pdMask);
%   pdShots(pdMask) = nan;

  ppeauShots = sum(pdShots, 1, 'omitnan'); % [1, nShots]
  % NOTE: taking max. does NOT correlate with power meter readings well!

  if nDim == 3
    ppeauShots = reshape(ppeauShots, [nx, ny]); % ppe in arbitrary unit
  end

end