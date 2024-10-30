%% script to get vessel maps with structural parameters like diameter and 
% angular alignment
clean_all;

data_base = 'C:\data_temp\woundExp';
woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day-3", "day4", "day7", "day10"];
timePt_list = ["day-3"];
nWound = length(woundID_list);
nTimePt = length(timePt_list);

%% Step 1: get GUI-compatible files for dermis structural maps

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)
        
        structure = load(fullfile(data_base, timePt, 'processed\structure', woundID));
        dermisMap = structure.structure_map_dermis;
        mapRaw = dermisMap;
        depthMap = ones(size(mapRaw));
        x = 1 : size(mapRaw,1); y = 1 : size(mapRaw,2);
        
        saveDir = fullfile(data_base, timePt, 'processed\structure\for_gui', woundID);
        save(saveDir, 'mapRaw', 'depthMap', 'x', 'y');
        
    end
    
end

%% Step 2: get binary masks and corresponding maps with structrual parameters

% TO BE CONTINUED in GUI: https://github.com/razanskylab/PostProGUI

woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)
        
        % read in processed vessel data
        vd = load(fullfile(data_base, timePt, 'processed\structure\for_gui\gui_export', strcat(woundID, "_vessels")));
        vList = vd.AVA.Data.vessel_list;
        nVessels = vd.AVA.nVessels;
        nVesselSegs = vd.AVA.nSegments;
        angleAlignSeg = zeros(nVesselSegs,1);
        
        % read in wound boundary and calculate geometric wound center
        woundBoundary = ReadImageJROI(fullfile(data_base, timePt, 'processed\forMasking', strcat(woundID, ".roi")));
        boundaryCoordinates = woundBoundary.mnCoordinates;
        woundGeoCenter = round(mean(boundaryCoordinates, 1));
        
        % calculate angular alignments w.r.t. wound center
        xCtr = woundGeoCenter(1);
        yCtr = woundGeoCenter(2);
        fun = @(x) cat(1, x);
        segCenter = cellfun(fun, {vList.centre}, 'UniformOutput', false);
        unitVectors = cellfun(fun, {vList.angles}, 'UniformOutput', false);
        nFills = 0;
        for iVessel = 1 : nVessels
            % calculate angle of segments should be -> centerAngle
            isegCenter = segCenter{iVessel};
            xDist = isegCenter(:,1)' - xCtr;
            yDist = isegCenter(:,2)' - yCtr;
            centerAngle = atan2d(xDist, yDist);
            centerAngle(centerAngle>90) = centerAngle(centerAngle>90) - 180; % only use +/- 90 deg
            centerAngle(centerAngle<-90) = centerAngle(centerAngle<-90) + 180; % only use +/- 90 deg
            % calculate angle of segments actually was -> segAngles
            iUnitVec = unitVectors{iVessel};
            segAngles = -atan2d(iUnitVec(:, 2), iUnitVec(:, 1))';
            % calculate difference between segAngles and centerAngle -> angleDiff
            angleDiff = segAngles(:) - centerAngle(:); 
            angleDiff(angleDiff>90) = 180 - angleDiff(angleDiff>90);
            angleDiff(angleDiff<-90) = 180 + angleDiff(angleDiff<-90); % only use +/- 90 deg
            angleDiff = abs(angleDiff);
            % 45 ==  mean(angleDiff) == random alignment
            % 0 == no diff -> full aligment
            % 90 == max diff, perpendicular aligned
            % convert 0-90 scale to 0-1 scale with 1 = full, 0 random and -1 missaling.
            angleAlign = (45-angleDiff)./45;
            nSeg = length(angleAlign);
            segStartIdx = nFills + 1;
            segEndIdx = segStartIdx + nSeg - 1;
            angleAlignSeg(segStartIdx:segEndIdx) = angleAlign;
%             averageAlignment(iVessel) = mean(angleAlign);

            nFills = nFills + nSeg;

        end
        
        % generate alignment maps and save
        vesselMask = vd.AVA.Data.bw';
        vesselSegCenter = round(vd.VesselData.segCenter);
        nPixels = size(vesselMask,1);
        alignMap = zeros(nPixels,nPixels);
        linearIdxSegCenter = sub2ind([nPixels,nPixels], vesselSegCenter(1,:), vesselSegCenter(2,:));
        alignMap(linearIdxSegCenter) = angleAlignSeg;
        
        % load wound masks
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        alignMap = fliplr(alignMap);
        alignMap = alignMap .* bandMasks.woundMask;
        
        % save
        saveDir = fullfile(data_base, timePt, 'processed\paramMaps\newAlignMaps', woundID);
        save(saveDir, 'alignMap');
        
    end
    
end

%% Step 3: save for further processing
% this section runs on the basis of exported vessel data from GUI

woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];
% timePt_list = ["day-3"];

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)
        
        % read in vessel parameters
        vd = load(fullfile(data_base, timePt, 'processed\structure\for_gui\gui_export', strcat(woundID, "_vessels")));
        vesselMask = vd.AVA.Data.bw';
        vesselSegCenter = round(vd.VesselData.segCenter);
        vesselCenter = round(vd.VesselData.vesCenter);
        vesselSegDiameter = vd.VesselData.segDiameters;
        vesselSegAlign = vd.VesselData.angleAlign;
        vesselTortuosity = vd.VesselData.turtosity;
        vesselAngleChange = vd.VesselData.angleChange;
        
        % generate vessel parameter maps
        nPixels = size(vesselMask,1);
        diameterMap = zeros(nPixels,nPixels);
        alignMap = zeros(nPixels,nPixels);
        tortuosityMap = zeros(nPixels,nPixels);
        angleChangeMap = zeros(nPixels,nPixels);
        linearIdxSegCenter = sub2ind([nPixels,nPixels], vesselSegCenter(1,:), vesselSegCenter(2,:));
        linearIdxVesCenter = sub2ind([nPixels,nPixels], vesselCenter(1,:), vesselCenter(2,:));
        diameterMap(linearIdxSegCenter) = vesselSegDiameter;
        alignMap(linearIdxSegCenter) = vesselSegAlign;
        tortuosityMap(linearIdxVesCenter) = vesselTortuosity;
        angleChangeMap(linearIdxVesCenter) = vesselAngleChange;
        
        % load wound masks
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        
        % make sure the parameter maps and wound masks align
        diameterMap = fliplr(diameterMap);
        alignMap = fliplr(alignMap);
        tortuosityMap = fliplr(tortuosityMap);
        angleChangeMap = fliplr(angleChangeMap);
        
        % mask out wound area in vessel parameter maps
        diameterMap = diameterMap .* bandMasks.woundMask;
        alignMap = alignMap .* bandMasks.woundMask;
        tortuosityMap = tortuosityMap .* bandMasks.woundMask;
        angleChangeMap = angleChangeMap .* bandMasks.woundMask;
        
        % save the vessel parameter maps for later statistics
        saveDir = fullfile(data_base, timePt, 'processed\paramMaps', woundID);
        save(saveDir, 'diameterMap', 'alignMap');

        saveDir = fullfile(data_base, timePt, 'processed\paramMaps\tortuosityMaps', woundID);
        save(saveDir, 'tortuosityMap', 'angleChangeMap');
        
    end
    
end

