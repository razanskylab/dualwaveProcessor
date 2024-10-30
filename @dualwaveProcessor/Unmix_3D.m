% File: Unmix_3D.m @ dualwaveProcessor
% Date: August 23, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Unmix_3D(dwp)

  fprintf('[dualwaveProcessor] Unmix in 3D, correcting for wavelength-depedent fluence attenuation...\n');

  % get necessary data and parameters
  usVol532 = dwp.usVol_532;
  usVol558 = dwp.usVol_dye;
  ppeCorrMap532 = dwp.ppeCorrectionMap_532;
  ppeCorrMap558 = dwp.ppeCorrectionMap_dye;

  mecMat = dwp.mecMat;
  vesselSigThreshold = dwp.vesselSigLB;
  maxNVessels = dwp.maxNVessels;

  skinCropDepths_mm_oa = dwp.tissueCropDepths_mm_oa;
  zVec_mm_532 = dwp.zVec_mm(1,:);
  zVec_mm_558 = dwp.zVec_mm(2,:);
  samplingDist_mm = dwp.samplingDist_mm;
  aboveSkin_mm = dwp.aboveTissue_mm;
  belowSkin_mm = dwp.belowTissue_mm;
  epidermisThickness_mm = dwp.epidermisThickness_mm;
  [~, nAlines_oa, nBscan_oa] = size(usVol532);

  % get already calculated scaling factors
  scaleFctIncident = dwp.scaleFctIncident; % incident PPE between-wavelength
  scaleFctReflectance = dwp.scaleFctReflectance; % surface reflectance

  % scattering coefficients in the skin
  miu_s_532 = get_reduced_scattering_coeff_skin(532);
  miu_s_558 = get_reduced_scattering_coeff_skin(558);

  % pre-calculate scale factor for attenuation through Epidermis
  if (epidermisThickness_mm == 0)
    warning('[dualwaveProcessor - Unmix3D] Assuming NO epidermis!');
    scaleFctEpidermis = 1;
  else
    miu_a_epi_532 = get_absorption_coeff_epidermis(532, dwp.Cm);
    miu_a_epi_558 = get_absorption_coeff_epidermis(558, dwp.Cm);
    phiEpidermis532 = estim_epidermis_attenuation(epidermisThickness_mm, miu_a_epi_532, miu_s_532);
    phiEpidermis558 = estim_epidermis_attenuation(epidermisThickness_mm, miu_a_epi_558, miu_s_558);
    scaleFctEpidermis = phiEpidermis532(end) / phiEpidermis558(end);
  end

  sO2_map = zeros(nAlines_oa, nBscan_oa, 'single');
  vesselDepth_map = zeros(nAlines_oa, nBscan_oa, 'single');
  vesselStructure_map = zeros(nAlines_oa, nBscan_oa, 'single');
  attenuation_map = zeros(nAlines_oa, nBscan_oa, 'single');

  % loop over each scanning position
  for iBscan = 1 : nBscan_oa
    if ~mod(iBscan, 100) % update progress every 100 B-scans
      fprintf('%d out of %d B-scans ...\n', iBscan, nBscan_oa);
    end
    for iAline = 1 : nAlines_oa
        
      currAline532 = usVol532(:,iAline,iBscan);
      currAline558 = usVol558(:,iAline,iBscan);
      currSkinDepth_mm = skinCropDepths_mm_oa(iAline, iBscan);
      
      % vessel detection
      [nVessel532, vesselIndices532] = detect_vessels_1D(currAline532, vesselSigThreshold, maxNVessels);
      [nVessel558, vesselIndices558] = detect_vessels_1D(currAline558, vesselSigThreshold, maxNVessels);
      
      if abs(vesselIndices558-vesselIndices532) > 1
        % the pair of A-lines at both wavelengths should NOT be depth shifted 
        continue;
      end
      
      if (nVessel532 == 0) || (nVessel558 == 0)
        continue;
      elseif (nVessel532 == 1) && (nVessel558 == 1)

        % extract vessel signal amplitude
        p0_532 = currAline532(vesselIndices532);
        p0_558 = currAline558(vesselIndices558);

        % correct for incident fluence
        p0_532 = p0_532 * ppeCorrMap532(iAline,iBscan);
        p0_558 = p0_558 * ppeCorrMap558(iAline,iBscan);
        p0_558 = p0_558 * scaleFctIncident;

        % correct for surface reflectance
        p0_558 = p0_558 * scaleFctReflectance;

        % correct for epidermis effective attenuation
        p0_558 = p0_558 * scaleFctEpidermis;

        % correct for dermis effective attenuation
        % the absorption in the dermis is dominated by blood, so in the
        % extra-vascular space, effective attenuation is dominated by scattering
        currVesselDepth_mm_532 = zVec_mm_532(vesselIndices532);
        currVesselDepth_mm_558 = zVec_mm_558(vesselIndices558);
        currDermisAttenuLength_mm_532 = currVesselDepth_mm_532 - (currSkinDepth_mm+epidermisThickness_mm);
        currDermisAttenuLength_mm_558 = currVesselDepth_mm_558 - (currSkinDepth_mm+epidermisThickness_mm);
        if (currDermisAttenuLength_mm_532<=samplingDist_mm) || (currDermisAttenuLength_mm_558<=samplingDist_mm)
            continue;
        end
        phiDermis532 = estim_dermis_attenuation(currDermisAttenuLength_mm_532, miu_s_532);
        phiDermis558 = estim_dermis_attenuation(currDermisAttenuLength_mm_558, miu_s_558);
        scaleFctDermis = phiDermis532(end) / phiDermis558(end);
        p0_558 = p0_558 * scaleFctDermis;

        p0_dw = [p0_532; p0_558];

        c_vec = mecMat \ p0_dw;
        curr_sO2 = c_vec(1) / sum(c_vec(:));

        sO2_map(iAline, iBscan) = curr_sO2;
        vesselDepth_map(iAline, iBscan) = currVesselDepth_mm_532;
        vesselStructure_map(iAline, iBscan) = p0_532;
        attenuation_map(iAline, iBscan) = scaleFctDermis;
          
      else
        error('Invalid number of detected vessels!');
      end
        
    end
  end

  if dwp.flagSNRMask
    maxVal = single(dwp.scanSett.sensitivityUs);
    amplitudeLowerThreshold = dwp.amplitudeLowerThresholdPct * maxVal;
    amplitudeUpperThreshold = dwp.amplitudeUpperThresholdPct * maxVal;
    snrMaskLower = (vesselStructure_map > amplitudeLowerThreshold);
    snrMaskUpper = (vesselStructure_map < amplitudeUpperThreshold);
    snrMask = snrMaskLower .* snrMaskUpper;
    sO2_map = sO2_map .* snrMask;
    vesselDepth_map = vesselDepth_map .* snrMask;
    vesselStructure_map = vesselStructure_map .* snrMask;
  end

  sO2_map(sO2_map==0) = nan;
  vesselDepth_map(vesselDepth_map==0) = nan;
  vesselStructure_map(vesselStructure_map==0) = nan;
  dwp.sO2_map = sO2_map;
  dwp.vesselDepth_map = vesselDepth_map;
  dwp.vesselStructure_map = vesselStructure_map;
  dwp.snrMask = snrMask;
  dwp.attenuation_map = attenuation_map;

  clear usVol532 usVol558;
  fprintf('done!\n');

end