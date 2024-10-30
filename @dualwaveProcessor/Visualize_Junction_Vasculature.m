% File: Visualize_Junction_Vasculature.m @ dualwaveProcessor
% Date: October 23, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function [structure_map_junction, sO2_map_junction] = Visualize_Junction_Vasculature(dwp)

  fprintf('[dualwaveProcessor] Visualizing vasculature at the dermis-hypodermis junction...');

  % getting necessary data and parameters
  vesselDepth_map = dwp.vesselDepth_map;
  sO2_map = dwp.sO2_map;
  vesselStructure_map = dwp.vesselStructure_map;
  skinCropDepths_mm_oa = dwp.tissueCropDepths_mm_oa;
  dermisEndThickness_mm = dwp.epidermisThickness_mm + dwp.dermisThickness_mm;
  junctionThickness_mm = dwp.junctionThickness_mm;
  [~, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);

  structure_map_junction = zeros(nAlines_oa, nBscan_oa);
  sO2_map_junction = zeros(nAlines_oa, nBscan_oa);

  % loop over each scanning position
  for iBscan = 1 : nBscan_oa
    if ~mod(iBscan, 100) % update progress every 100 B-scans
      fprintf('%d out of %d B-scans ...\n', iBscan, nBscan_oa);
    end
    for iAline = 1 : nAlines_oa

      currVesselDepth_mm = vesselDepth_map(iAline,iBscan);
      currSkinDepth_mm = skinCropDepths_mm_oa(iAline,iBscan);
      currJunctionDepth_mm = currSkinDepth_mm + dermisEndThickness_mm;
      curr_sO2 = sO2_map(iAline,iBscan);
      curr_structure = vesselStructure_map(iAline,iBscan);
      if (isnan(currVesselDepth_mm) || (currVesselDepth_mm<currSkinDepth_mm))
        continue;
      end

      if ((currVesselDepth_mm>(currJunctionDepth_mm-junctionThickness_mm/2))&&(currVesselDepth_mm<(currJunctionDepth_mm+junctionThickness_mm/2)))
        % this vessel belongs to the junction slice
        structure_map_junction(iAline,iBscan) = curr_structure;
        sO2_map_junction(iAline,iBscan) = curr_sO2;
      else
        continue;
      end

    end
  end

end