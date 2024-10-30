% File: Render_sO2_3D.m @ dualwaveProcessor
% Date: August 25, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function Render_sO2_3D(dwp)

  fprintf('[dualwaveProcessor] Rendering 3D structural and sO2 images...');

  % getting necessary data and parameters
  vesselDepth_map = dwp.vesselDepth_map;
  sO2_map = dwp.sO2_map;

  zVec_mm_532 = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
  zVec_mm_532 = zVec_mm_532(dwp.zCropIdx_us(1,1):dwp.zCropIdx_us(1,2));

  [nt, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);
  sO2_vol = zeros(nt, nAlines_oa, nBscan_oa, 'single');

  % loop over every scanning position
  for iBscan = 1 : nBscan_oa
    if ~mod(iBscan, 100) % update progress every 100 B-scans
      fprintf('%d out of %d B-scans ...\n', iBscan, nBscan_oa);
    end
      for iAline = 1 : nAlines_oa

        curr_sO2 = sO2_map(iAline, iBscan);
        currVesselDepth_mm = vesselDepth_map(iAline, iBscan);
        if isnan(curr_sO2)
          continue;
        end

        % assign sO2 value to vessel voxels along depth
        currAline = dwp.usVol_532(:,iAline,iBscan);
        currAline(currAline<0) = 0;
        [~,currVesselIdx] = min(abs(zVec_mm_532-currVesselDepth_mm));
        zeroIndices = find(~currAline);
        [~,closestZeroIndices] = mink(abs(currVesselIdx-zeroIndices), 2);
        currVesselIndices = zeroIndices(closestZeroIndices);
        currVesselIndices = sort(currVesselIndices, 'ascend');
        sO2_vol(currVesselIndices(1):currVesselIndices(2), iAline, iBscan) = curr_sO2;
          
      end
  end

  % assign to property
  dwp.sO2_vol = sO2_vol;

  fprintf('done!\n');

end