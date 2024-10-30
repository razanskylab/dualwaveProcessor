% File: Extract_Tissue_Data.m @ dualwaveProcessor
% Date: August 25, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Extract_Tissue_Data(dwp)

  fprintf('[dualwaveProcessor] Extracting tissue data based on estimated tissue surface...');

  % getting necessary data and parameters
  skinCropDepths_mm_oa = dwp.tissueCropDepths_mm_oa;
  zVec_mm_532 = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
  zVec_mm_532 = zVec_mm_532(dwp.zCropIdx_us(1,1):dwp.zCropIdx_us(1,2));
  zVec_mm_558 = dwp.scanSett.zVec(2,:) + dwp.scanSett.fd*1e3;
  zVec_mm_558 = zVec_mm_558(dwp.zCropIdx_us(2,1):dwp.zCropIdx_us(2,2));
  [nt, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);
  samplingDist_mm = dwp.samplingDist_mm;
  aboveSkin_mm = dwp.aboveTissue_mm;
  belowSkin_mm = dwp.belowTissue_mm;
  % allocate memory for tissue data
  nPixel_abveSkin = aboveSkin_mm / samplingDist_mm;
  nPixel_belowSkin = belowSkin_mm / samplingDist_mm;
  nPixel_allSkin = nPixel_abveSkin + nPixel_belowSkin;
  usVol532 = zeros(nPixel_allSkin, nAlines_oa, nBscan_oa);
  usVol558 = zeros(nPixel_allSkin, nAlines_oa, nBscan_oa);
  sO2Vol = zeros(nPixel_allSkin, nAlines_oa, nBscan_oa);

  % loop over every scanning position
  for iBscan = 1 : nBscan_oa
    if ~mod(iBscan, 100) % update progress every 100 B-scans
      fprintf('%d out of %d B-scans ...\n', iBscan, nBscan_oa);
    end
      for iAline = 1 : nAlines_oa
          
        currSkinDepth_mm = skinCropDepths_mm_oa(iAline, iBscan);
        [~,currCropIdx532] = min(abs(zVec_mm_532 - currSkinDepth_mm));
        [~,currCropIdx558] = min(abs(zVec_mm_558 - currSkinDepth_mm));
        startIdx532 = currCropIdx532 - nPixel_abveSkin + 1;
        endIdx532 = currCropIdx532 + nPixel_belowSkin;
        startIdx558 = currCropIdx558 - nPixel_abveSkin + 1;
        endIdx558 = currCropIdx558 + nPixel_belowSkin;
        usVol532(:,iAline,iBscan) = dwp.usVol_532(startIdx532:endIdx532,iAline,iBscan);
        usVol558(:,iAline,iBscan) = dwp.usVol_dye(startIdx558:endIdx558,iAline,iBscan);
        sO2Vol(:,iAline,iBscan) = dwp.sO2_vol(startIdx532:endIdx532,iAline,iBscan);
          
      end
  end

  % assign to property
  dwp.usVolTissue_532 = usVol532;
  dwp.usVolTissue_dye = usVol558;
  dwp.sO2VolTissue = sO2Vol;
  clear usVol532 usVol558 sO2Vol;

  fprintf('done!\n');

end