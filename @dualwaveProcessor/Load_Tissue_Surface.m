% File: Load_Tissue_Surface.m @ dualwaveProcessor
% Date: August 23, 2023
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

% NOTE: now only skin surface is implemented

function Load_Tissue_Surface(dwp)

  fprintf('[dualwaveProcessor] Loading estimated tissue surface...');

  % load estimated tissue surface from US dataset
  currDataName = dwp.data_name;
  mouseInfo = currDataName(1:8);
  skinCropData = matfile(fullfile(dwp.data_dir, 'skinCropDepths', mouseInfo));
  skinCropDepths_mm_us = skinCropData.skinCropDepths_mm;
  [nAlines_us, nBscan_us] = size(skinCropDepths_mm_us);
  % assign to property
  dwp.tissueCropDepths_mm_us = skinCropDepths_mm_us;
  
  % get US scanning grid
  stepSize_mm_us = dwp.stepSizeWound_mm_us;
  xVec_mm_us = (0:nAlines_us-1) .* stepSize_mm_us;
  yVec_mm_us = (0:nBscan_us-1) .* stepSize_mm_us;
  usGrid{1} = xVec_mm_us;
  usGrid{2} = yVec_mm_us;

  % get OA scanning grid
  [~, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);
  stepSize_mm_oa = dwp.stepSizeWound_mm_oa;
  xVec_mm_oa = (0:nAlines_oa-1) .* stepSize_mm_oa;
  yVec_mm_oa = (0:nBscan_oa-1) .* stepSize_mm_oa;
  oaGrid{1} = xVec_mm_oa;
  oaGrid{2} = yVec_mm_oa;

  interpF = griddedInterpolant(usGrid, skinCropDepths_mm_us, 'linear', 'linear');
  skinCropDepths_mm_oa = interpF(oaGrid);

  % [~, nAlinesZoomin_oa, nBscanZoomin_oa] = size(dwp.usVol_532);
  % stepSize_mm_zoomin = dwp.stepSize_mm_zoomin;
  % FOVSize_mm_zoomin = (nAlinesZoomin_oa-1) * stepSize_mm_zoomin;

  % % If the current dataset is a zoom-in scan over a sub-FOV, select the 
  % % corresponding area in skin surface and interpolate to even finer zoom-in scanning grid
  % if dwp.flagZoominScan
  %   zoominFOVCenter_mm = dwp.scanSett.center;
  %   dwp.zoominFOVCenter = zoominFOVCenter_mm;
  %   fullFOVCenter_mm = dwp.fullFOVCenter;
  %   deltaCenter_mm = zoominFOVCenter_mm - fullFOVCenter_mm;
  %   deltaCenter_idx = deltaCenter_mm ./ stepSize_mm_oa;
  %   fullFOVCenter_idx = round([nAlinesFull_oa, nBscanFull_oa]./2);
  %   zoominFOVCenter_idx = fullFOVCenter_idx + flip(deltaCenter_idx); % [x,y] -> [nAlines, nBscans]
  %   % prepare interpolation grid
  %   FOVSize_idx_fullFOVStep = round(FOVSize_mm_zoomin / stepSize_mm_oa);
  %   FOVSize_idx_zoominFOVStep = round(FOVSize_mm_zoomin / stepSize_mm_zoomin);
  %   fullGrid{1} = (0:FOVSize_idx_fullFOVStep) .* stepSize_mm_oa;
  %   fullGrid{2} = (0:FOVSize_idx_fullFOVStep) .* stepSize_mm_oa;
  %   zoominGrid{1} = (0:FOVSize_idx_zoominFOVStep) .* stepSize_mm_zoomin;
  %   zoominGrid{2} = (0:FOVSize_idx_zoominFOVStep) .* stepSize_mm_zoomin;
  %   skinCropDepths_mm_fullFOVStep = skinCropDepths_mm_oa(...
  %                                     (zoominFOVCenter_idx(1)-(FOVSize_idx_fullFOVStep/2)) : (zoominFOVCenter_idx(1)+(FOVSize_idx_fullFOVStep/2)), ...
  %                                     (zoominFOVCenter_idx(2)-(FOVSize_idx_fullFOVStep/2)) : (zoominFOVCenter_idx(2)+(FOVSize_idx_fullFOVStep/2)));
  %   interpF_zoomin = griddedInterpolant(fullGrid, skinCropDepths_mm_fullFOVStep, 'linear', 'linear');
  %   skinCropDepths_mm_zoominFOVStep = interpF_zoomin(zoominGrid);
  %   % check dimensions of skin surface profile and actual dataset match
  %   if size(skinCropDepths_mm_zoominFOVStep,1) ~= nAlinesZoomin_oa
  %     error('Number of A-lines does NOT match between tissue surface and actual Zoom-in dataset!');
  %   end
  %   if size(skinCropDepths_mm_zoominFOVStep,2) ~= nBscanZoomin_oa
  %     warning('Number of B-scans does NOT match between tissue surface and actual Zoom-in dataset!');
  %     warning('Getting grid of the last B-scan...');
  %     skinCropDepths_mm_zoominFOVStep = skinCropDepths_mm_zoominFOVStep(:,1:end-1);
  %     if size(skinCropDepths_mm_zoominFOVStep,2) ~= nBscanZoomin_oa
  %       error('Number of B-scans does NOT match between tissue surface and actual Zoom-in dataset after correction!');
  %     end
  %   end
  %   % assign to property
  %   dwp.tissueCropDepths_mm_oa = skinCropDepths_mm_zoominFOVStep;
  % else
  %   % assign to property
  %   dwp.tissueCropDepths_mm_oa = skinCropDepths_mm_oa;
  % end

  dwp.tissueCropDepths_mm_oa = skinCropDepths_mm_oa;

  clear skinCropDepths_mm_us skinCropDepths_mm_oa;
  clear skinCropDepths_mm_fullFOVStep skinCropDepths_mm_zoominFOVStep;

  fprintf('done!\n');

end