% File: Assign_To_Tissue_Layer.m @ dualwaveProcessor
% Date: October 10, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

function [mapDermis, mapHypodermis] = Assign_To_Tissue_Layer(dwp)

  fprintf('[dualwaveProcessor] Assigning vessels to their respective tissue layer...');

  % getting necessary data and parameters
  vesselDepth_map = dwp.vesselDepth_map;
  sO2_map = dwp.sO2_map;
  vesselStructure_map = dwp.vesselStructure_map;
  skinCropDepths_mm_oa = dwp.tissueCropDepths_mm_oa;
  dermisEndThickness_mm = dwp.epidermisThickness_mm + dwp.dermisThickness_mm;
  [~, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);

  mapDermis = zeros(nAlines_oa, nBscan_oa);
  mapHypodermis = zeros(nAlines_oa, nBscan_oa);

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

      % assign sO2 info
      if strcmp(dwp.assignToTissueLayerInfo, 'sO2')
        if (currVesselDepth_mm < currJunctionDepth_mm)
          mapDermis(iAline,iBscan) = curr_sO2;
        else
          mapHypodermis(iAline,iBscan) = curr_sO2;
        end
      end

      % assign structure info
      if strcmp(dwp.assignToTissueLayerInfo, 'structure')
        if (currVesselDepth_mm < currJunctionDepth_mm)
          mapDermis(iAline,iBscan) = curr_structure;
        else
          mapHypodermis(iAline,iBscan) = curr_structure;
        end
      end

    end
  end

end