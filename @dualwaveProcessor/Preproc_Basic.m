% File: Preproc_Basic.m @ dualwaveProcessor
% Date: July 28th, 2022
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Preproc_Basic(dwp)

  fprintf('[dualwaveProcessor] Basic pre-processing including jitter, flipping, cropping, coverting data type and averaging...');

  mode = dwp.scanSett.mode;

  [~, nWavelength, nAverage, ~, ~] = size(dwp.RawDataUs);
  if nAverage > 1
    error('Averaging function NOT ready yet!');
  else
    usData = squeeze(dwp.RawDataUs); % [ntUs, nWavelength, nx, ny]
    pdData = squeeze(dwp.RawDataPd); % [ntPd, nWavelength, nx, ny]
  end

  if strcmp(mode, 'dualwave')
    assert(nWavelength == 2);
    % split wavelength
    usVol_532 = squeeze(usData(:,1,:,:)); % [ntUs, nx, ny]
    usVol_dye = squeeze(usData(:,2,:,:)); % [ntUs, nx, ny]
    pdVol_532 = squeeze(pdData(:,1,:,:)); % [ntPd, nx, ny]
    pdVol_dye = squeeze(pdData(:,2,:,:)); % [ntPd, nx, ny]
    % correct for sinusoidal scan pattern
    usVol_532(:,:,2:2:end) = flip(usVol_532(:,:,2:2:end), 2);
    usVol_dye(:,:,2:2:end) = flip(usVol_dye(:,:,2:2:end), 2);
    pdVol_532(:,:,2:2:end) = flip(pdVol_532(:,:,2:2:end), 2);
    pdVol_dye(:,:,2:2:end) = flip(pdVol_dye(:,:,2:2:end), 2);
    % jitter correction
    [pdVol_532, usVol_532, pdMaxIdx_532] = correct_laser_jitter(pdVol_532, usVol_532);
    [pdVol_dye, usVol_dye, pdMaxIdx_dye] = correct_laser_jitter(pdVol_dye, usVol_dye);
    % update PD initial peaks for accurate calculation of depth vector
    dwp.scanSett.pdInitPeak_532 = pdMaxIdx_532;
    dwp.scanSett.pdInitPeak_dye = pdMaxIdx_dye;
    % covert DAQ count to voltage
    usVol_532 = counts_to_voltage(usVol_532, dwp.scanSett.sensitivityUs, dwp.daqResolution);
    usVol_dye = counts_to_voltage(usVol_dye, dwp.scanSett.sensitivityUs, dwp.daqResolution);
    pdVol_532 = counts_to_voltage(pdVol_532, dwp.scanSett.sensitivityPd, dwp.daqResolution);
    pdVol_dye = counts_to_voltage(pdVol_dye, dwp.scanSett.sensitivityPd, dwp.daqResolution);
    % crop along depth
    zCropIdx_us = dwp.zCropIdx_us;
    usVol_532 = usVol_532(zCropIdx_us(1,1):zCropIdx_us(1,2), :, :);
    usVol_dye = usVol_dye(zCropIdx_us(2,1):zCropIdx_us(2,2), :, :);
    % assign to property
    dwp.usVol_532 = usVol_532;
    dwp.usVol_dye = usVol_dye;
    dwp.pdVol_532 = pdVol_532;
    dwp.pdVol_dye = pdVol_dye;

  elseif strcmp(mode, 'onda532')
    assert(nWavelength == 1);
    % jitter correction
    [pdData, usData, pdMaxIdx_532] = correct_laser_jitter(pdData, usData);
    % update PD initial peaks for accurate calculation of depth vector
    dwp.scanSett.pdInitPeak_532 = pdMaxIdx_532;
    % correct for sinusoidal scan pattern
    usData(:,:,2:2:end) = flip(usData(:,:,2:2:end), 2);
    pdData(:,:,2:2:end) = flip(pdData(:,:,2:2:end), 2);
    % covert DAQ count to voltage
    usData = counts_to_voltage(usData, dwp.scanSett.sensitivityUs, dwp.daqResolution);
    pdData = counts_to_voltage(pdData, dwp.scanSett.sensitivityPd, dwp.daqResolution);
    % crop along depth
    zCropIdx_us = dwp.zCropIdx_us;
    usData = usData(zCropIdx_us(1,1):zCropIdx_us(1,2), :, :);
    % assign to property
    dwp.usVol_532 = usData;
    dwp.pdVol_532 = pdData;

  elseif strcmp(mode, 'dye')
    assert(nWavelength == 1);
    % jitter correction
    [pdData, usData, pdMaxIdx_dye] = correct_laser_jitter(pdData, usData);
    % update PD initial peaks for accurate calculation of depth vector
    dwp.scanSett.pdInitPeak_dye = pdMaxIdx_dye;
    % correct for sinusoidal scan pattern
    usData(:,:,2:2:end) = flip(usData(:,:,2:2:end), 2);
    pdData(:,:,2:2:end) = flip(pdData(:,:,2:2:end), 2);
    % covert DAQ count to voltage
    usData = counts_to_voltage(usData, dwp.scanSett.sensitivityUs, dwp.daqResolution);
    pdData = counts_to_voltage(pdData, dwp.scanSett.sensitivityPd, dwp.daqResolution);
    % crop along depth
    zCropIdx_us = dwp.zCropIdx_us;
    usData = usData(zCropIdx_us(2,1):zCropIdx_us(2,2), :, :);

    % assign to property
    dwp.usVol_dye = usData;
    dwp.pdVol_dye = pdData;

  elseif strcmp(mode, 'us')
    assert(nWavelength == 1);
    % correct for sinusoidal scan pattern
    usData(:,:,2:2:end) = flip(usData(:,:,2:2:end), 2);
    % covert DAQ count to voltage
    usData = counts_to_voltage(usData, dwp.scanSett.sensitivityUs, dwp.daqResolution);
    % crop along depth
    zCropIdx_us = dwp.zCropIdx_us;
    usData = usData(zCropIdx_us(1,1):zCropIdx_us(1,2), :, :);
    % assign to property
    dwp.usVol_pe = usData;

  else
    error('[dualwaveProcessor] invalid scan mode!');
  end

  % free up space occupied by RawDataUs and RawDataPd
  dwp.RawDataPd = [];
  dwp.RawDataUs = [];

  fprintf('done!\n');

end