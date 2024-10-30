% File: Unmix_2D.m @ dualwaveProcessor
% Date: July 28th, 2022
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

% Change log on Feb. 19, 2023: disabled spectral slope masking
% Change log on August 25, 2023: enable 2D unmixing based on max. amplitude spectrum with the updated properties and methods

function Unmix_2D(dwp)

  fprintf('[dualwaveProcessor] Computing sO2 based on 2D processing of volumes...\n');

  mecMat = dwp.mecMat;
  scaleFctIncident = dwp.scaleFctIncident;

  [~, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);
  sO2_map = zeros(nAlines_oa, nBscan_oa, 'single');

  % estimate sO2 on an Aline-by-Aline basis
  for iBscan = 1 : nBscan_oa
    if ~mod(iBscan, 100) % update progress every 100 B-scans
      fprintf('%d out of %d B-scans ...\n', iBscan, nBscan_oa);
    end
    for iAline = 1 : nAlines_oa

      currAline_532 = dwp.usVol_532(:, iAline, iBscan);
      currAline_dye = dwp.usVol_dye(:, iAline, iBscan);

      % scale each pair of dual-wave OA signals here
      ppeScaleFct_532 = dwp.ppeCorrectionMap_532(iAline, iBscan);
      ppeScaleFct_dye = dwp.ppeCorrectionMap_dye(iAline, iBscan);
      if (isnan(ppeScaleFct_dye) || isnan(ppeScaleFct_532))
        continue; % jump to the next position directly
      else
        % Normalize to corresponding median PPE level at each wavelength
        currAline_532 = currAline_532 .* ppeScaleFct_532;
        currAline_dye = currAline_dye .* ppeScaleFct_dye;
        % Remove between-wavelength PPE difference
        currAline_dye = currAline_dye .* scaleFctIncident;
      end

      % calculate peak-to-peak for onda532
      [maxAmp_532, ~] = max(currAline_532);
      [minAmp_532, ~] = min(currAline_532);
      p0_532 = maxAmp_532 - minAmp_532;
      % calculate peak-to-peak for dye
      [maxAmp_dye, ~] = max(currAline_dye);
      [minAmp_dye, ~] = min(currAline_dye);
      p0_dye = maxAmp_dye - minAmp_dye;

      % use only max. amplitude
      p0_532 = maxAmp_532;
      p0_dye = maxAmp_dye;

      p0_dw = [p0_532; p0_dye];
      c_vec = mecMat \ p0_dw;
      curr_sO2 = c_vec(1) / sum(c_vec(:));

      sO2_map(iAline, iBscan) = curr_sO2;

    end
  end

  if dwp.flagSNRMask
    mip = squeeze(max(dwp.usVol_dye, [], 1));
    % maxVal = max(abs(mip(:)));
    maxVal = single(dwp.scanSett.sensitivityUs);
    amplitudeLowerThreshold = dwp.amplitudeLowerThresholdPct * maxVal;
    amplitudeUpperThreshold = dwp.amplitudeUpperThresholdPct * maxVal;
    snrMaskLower = (mip > amplitudeLowerThreshold);
    snrMaskUpper = (mip < amplitudeUpperThreshold);
    snrMask = snrMaskLower .* snrMaskUpper;
    sO2_map = sO2_map .* snrMask;
  end

  sO2_map(sO2_map == 0) = nan;
  dwp.sO2_map = sO2_map;
  dwp.snrMask = snrMask;
  fprintf('done!\n');

end