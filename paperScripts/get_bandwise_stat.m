%% get band-based vessel statistics, analyzed vascular parameters include:
% sO2, diameter, angular alignment and vessel density
clean_all;

data_base = 'C:\data_temp\woundExp';
woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];

nWound = length(woundID_list);
nTimePt = length(timePt_list);
nBand = 5;
medfiltSize = 3;

sO2_stat = zeros(nWound, nTimePt, nBand);
sO2_hetero_stat = zeros(nWound, nTimePt, nBand);
density_stat = zeros(nWound, nTimePt, nBand);
diameter_stat = zeros(nWound, nTimePt, nBand);
align_stat = zeros(nWound, nTimePt, nBand);
tortuosity_stat = zeros(nWound, nTimePt, nBand);
angleChange_stat = zeros(nWound, nTimePt, nBand);

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)

        structure = load(fullfile(data_base, timePt, 'processed\structure', woundID));
        sO2 = load(fullfile(data_base, timePt, 'processed\sO2', woundID));
        structuralParam = load(fullfile(data_base, timePt, 'processed\paramMaps', woundID));
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        newAlignment = load(fullfile(data_base, timePt, 'processed\paramMaps\newAlignMaps', woundID));
        tortuosityParam = load(fullfile(data_base, timePt, 'processed\paramMaps\tortuosityMaps', woundID));
        
        % structural information
        structure_map = structure.structure_map_dermis';
        structure_map = medfilt2(structure_map.*bandMasks.woundMask,[medfiltSize,medfiltSize]);
        diameter_map = structuralParam.diameterMap;
%         align_map = structuralParam.alignMap;
        align_map = newAlignment.alignMap;
        tortuosity_map = tortuosityParam.tortuosityMap;
        angleChange_map = tortuosityParam.angleChangeMap;

        % oxygenation information
        sO2_map = sO2.sO2_map_dermis';
        sO2_map = medfilt2(sO2_map.*bandMasks.woundMask,[medfiltSize,medfiltSize]);
        % remove nonsense sO2 values
        sO2_map(sO2_map>1) = 0;
        sO2_map(sO2_map<0) = 0;
        
        % band 1
        band1Mask = bandMasks.band1Mask;
        structure_band1 = structure_map .* band1Mask;       
        sO2_band1 = sO2_map .* band1Mask;
        density_band1 = nnz(structure_band1(:)) / nnz(band1Mask(:));
        [~,~,sO2Vals_band1] = find(sO2_band1(:));
        median_sO2_band1 = median(sO2Vals_band1(:), 'omitnan');
        sO2_hetero_band1 = std(sO2Vals_band1(:), 'omitnan');
        diameter_band1 = diameter_map .* band1Mask;
        [~,~,diameterVals_band1] = find(diameter_band1(:));
        median_diameter_band1 = median(diameterVals_band1(:), 'omitnan');
        align_band1 = align_map .* band1Mask;
        [~,~,alignVals_band1] = find(align_band1(:));
        median_align_band1 = median(alignVals_band1(:), 'omitnan');
        tortuosity_band1 = tortuosity_map .* band1Mask;
        [~,~,tortuosityVals_band1] = find(tortuosity_band1(:));
        median_tortuosity_band1 = median(tortuosityVals_band1(:), 'omitnan');
        angleChange_band1 = angleChange_map .* band1Mask;
        [~,~,angleChangeVals_band1] = find(angleChange_band1(:));
        median_angleChange_band1 = median(angleChangeVals_band1(:), 'omitnan');
        
        % band 2
        band2Mask = bandMasks.band2Mask;
        structure_band2 = structure_map .* band2Mask;
        sO2_band2 = sO2_map .* band2Mask;
        density_band2 = nnz(structure_band2(:)) / nnz(band2Mask(:));
        [~,~,sO2Vals_band2] = find(sO2_band2(:));
        median_sO2_band2 = median(sO2Vals_band2(:), 'omitnan');
        sO2_hetero_band2 = std(sO2Vals_band2(:), 'omitnan');
        diameter_band2 = diameter_map .* band2Mask;
        [~,~,diameterVals_band2] = find(diameter_band2(:));
        median_diameter_band2 = median(diameterVals_band2(:), 'omitnan');
        align_band2 = align_map .* band2Mask;
        [~,~,alignVals_band2] = find(align_band2(:));
        median_align_band2 = median(alignVals_band2(:), 'omitnan');
        tortuosity_band2 = tortuosity_map .* band2Mask;
        [~,~,tortuosityVals_band2] = find(tortuosity_band2(:));
        median_tortuosity_band2 = median(tortuosityVals_band2(:), 'omitnan');
        angleChange_band2 = angleChange_map .* band2Mask;
        [~,~,angleChangeVals_band2] = find(angleChange_band2(:));
        median_angleChange_band2 = median(angleChangeVals_band2(:), 'omitnan');
        
        % band 3
        band3Mask = bandMasks.band3Mask;
        structure_band3 = structure_map .* band3Mask;
        sO2_band3 = sO2_map .* band3Mask;
        density_band3 = nnz(structure_band3(:)) / nnz(band3Mask(:));
        [~,~,sO2Vals_band3] = find(sO2_band3(:));
        median_sO2_band3 = median(sO2Vals_band3(:), 'omitnan');
        sO2_hetero_band3 = std(sO2Vals_band3(:), 'omitnan');
        diameter_band3 = diameter_map .* band3Mask;
        [~,~,diameterVals_band3] = find(diameter_band3(:));
        median_diameter_band3 = median(diameterVals_band3(:), 'omitnan');
        align_band3 = align_map .* band3Mask;
        [~,~,alignVals_band3] = find(align_band3(:));
        median_align_band3 = median(alignVals_band3(:), 'omitnan');
        tortuosity_band3 = tortuosity_map .* band3Mask;
        [~,~,tortuosityVals_band3] = find(tortuosity_band3(:));
        median_tortuosity_band3 = median(tortuosityVals_band3(:), 'omitnan');
        angleChange_band3 = angleChange_map .* band3Mask;
        [~,~,angleChangeVals_band3] = find(angleChange_band3(:));
        median_angleChange_band3 = median(angleChangeVals_band3(:), 'omitnan');
        
        % band 4
        band4Mask = bandMasks.band4Mask;
        structure_band4 = structure_map .* band4Mask;
        sO2_band4 = sO2_map .* band4Mask;
        density_band4 = nnz(structure_band4(:)) / nnz(band4Mask(:));
        [~,~,sO2Vals_band4] = find(sO2_band4(:));
        median_sO2_band4 = median(sO2Vals_band4(:), 'omitnan');
        sO2_hetero_band4 = std(sO2Vals_band4(:), 'omitnan');
        diameter_band4 = diameter_map .* band4Mask;
        [~,~,diameterVals_band4] = find(diameter_band4(:));
        median_diameter_band4 = median(diameterVals_band4(:), 'omitnan');
        align_band4 = align_map .* band4Mask;
        [~,~,alignVals_band4] = find(align_band4(:));
        median_align_band4 = median(alignVals_band4(:), 'omitnan');
        tortuosity_band4 = tortuosity_map .* band4Mask;
        [~,~,tortuosityVals_band4] = find(tortuosity_band4(:));
        median_tortuosity_band4 = median(tortuosityVals_band4(:), 'omitnan');
        angleChange_band4 = angleChange_map .* band4Mask;
        [~,~,angleChangeVals_band4] = find(angleChange_band4(:));
        median_angleChange_band4 = median(angleChangeVals_band4(:), 'omitnan');
        
        % band 5
        band5Mask = bandMasks.band5Mask;
        structure_band5 = structure_map .* band5Mask;
        sO2_band5 = sO2_map .* band5Mask;
        density_band5 = nnz(structure_band5(:)) / nnz(band5Mask(:));
        [~,~,sO2Vals_band5] = find(sO2_band5(:));
        median_sO2_band5 = median(sO2Vals_band5(:), 'omitnan');
        sO2_hetero_band5 = std(sO2Vals_band5(:), 'omitnan');
        diameter_band5 = diameter_map .* band5Mask;
        [~,~,diameterVals_band5] = find(diameter_band5(:));
        median_diameter_band5 = median(diameterVals_band5(:), 'omitnan');
        align_band5 = align_map .* band5Mask;
        [~,~,alignVals_band5] = find(align_band5(:));
        median_align_band5 = median(alignVals_band5(:), 'omitnan');
        tortuosity_band5 = tortuosity_map .* band5Mask;
        [~,~,tortuosityVals_band5] = find(tortuosity_band5(:));
        median_tortuosity_band5 = median(tortuosityVals_band5(:), 'omitnan');
        angleChange_band5 = angleChange_map .* band5Mask;
        [~,~,angleChangeVals_band5] = find(angleChange_band5(:));
        median_angleChange_band5 = median(angleChangeVals_band5(:), 'omitnan');
        
        % assign to sO2 values
        sO2_stat(iWound, iTimePt, 1) = median_sO2_band1;
        sO2_stat(iWound, iTimePt, 2) = median_sO2_band2;
        sO2_stat(iWound, iTimePt, 3) = median_sO2_band3;
        sO2_stat(iWound, iTimePt, 4) = median_sO2_band4;
        sO2_stat(iWound, iTimePt, 5) = median_sO2_band5;
        
        % assign to sO2 heterogeneity values
        sO2_hetero_stat(iWound, iTimePt, 1) = sO2_hetero_band1;
        sO2_hetero_stat(iWound, iTimePt, 2) = sO2_hetero_band2;
        sO2_hetero_stat(iWound, iTimePt, 3) = sO2_hetero_band3;
        sO2_hetero_stat(iWound, iTimePt, 4) = sO2_hetero_band4;
        sO2_hetero_stat(iWound, iTimePt, 5) = sO2_hetero_band5;
        
        % assign to density values
        density_stat(iWound, iTimePt, 1) = density_band1;
        density_stat(iWound, iTimePt, 2) = density_band2;
        density_stat(iWound, iTimePt, 3) = density_band3;
        density_stat(iWound, iTimePt, 4) = density_band4;
        density_stat(iWound, iTimePt, 5) = density_band5;
        
        % assign to diameter values
        diameter_stat(iWound, iTimePt, 1) = median_diameter_band1;
        diameter_stat(iWound, iTimePt, 2) = median_diameter_band2;
        diameter_stat(iWound, iTimePt, 3) = median_diameter_band3;
        diameter_stat(iWound, iTimePt, 4) = median_diameter_band4;
        diameter_stat(iWound, iTimePt, 5) = median_diameter_band5;
        
        % assign to angular alignment values
        align_stat(iWound, iTimePt, 1) = median_align_band1;
        align_stat(iWound, iTimePt, 2) = median_align_band2;
        align_stat(iWound, iTimePt, 3) = median_align_band3;
        align_stat(iWound, iTimePt, 4) = median_align_band4;
        align_stat(iWound, iTimePt, 5) = median_align_band5;
        
        % assign to tortuosity values
        tortuosity_stat(iWound, iTimePt, 1) = median_tortuosity_band1;
        tortuosity_stat(iWound, iTimePt, 2) = median_tortuosity_band2;
        tortuosity_stat(iWound, iTimePt, 3) = median_tortuosity_band3;
        tortuosity_stat(iWound, iTimePt, 4) = median_tortuosity_band4;
        tortuosity_stat(iWound, iTimePt, 5) = median_tortuosity_band5;
        
        % assign to angleChange values
        angleChange_stat(iWound, iTimePt, 1) = median_angleChange_band1;
        angleChange_stat(iWound, iTimePt, 2) = median_angleChange_band2;
        angleChange_stat(iWound, iTimePt, 3) = median_angleChange_band3;
        angleChange_stat(iWound, iTimePt, 4) = median_angleChange_band4;
        angleChange_stat(iWound, iTimePt, 5) = median_angleChange_band5;
        
    end
end

% convert to proper units
sO2_stat = sO2_stat .* 100; % [%]
sO2_hetero_stat = sO2_hetero_stat .* 100; % [%]
diameter_stat = diameter_stat .* 5; % [um]
% align_stat(align_stat<0) = nan; % get rid of outlier median angular alignment values

%% show band-wise sO2, diameter and angular alignment heatmaps
% NOTE: works for one selected wound

% sO2
sO2Heatmap = zeros(size(band1Mask));
sO2Heatmap = sO2Heatmap + band1Mask.*median_sO2_band1;
sO2Heatmap = sO2Heatmap + band2Mask.*median_sO2_band2;
sO2Heatmap = sO2Heatmap + band3Mask.*median_sO2_band3;
sO2Heatmap = sO2Heatmap + band4Mask.*median_sO2_band4;
sO2Heatmap = sO2Heatmap + band5Mask.*median_sO2_band5;
figure('Name', 'sO2 heatmap'); set(gcf, 'color', 'w');
imagesc(sO2Heatmap); axis image; colormap inferno;
caxis([30 80]); xticks([]); yticks([]);
h = colorbar; h.Ticks = [];

% diameter
diamHeatmap = zeros(size(band1Mask));
diamHeatmap = diamHeatmap + band1Mask.*median_diameter_band1;
diamHeatmap = diamHeatmap + band2Mask.*median_diameter_band2;
diamHeatmap = diamHeatmap + band3Mask.*median_diameter_band3;
diamHeatmap = diamHeatmap + band4Mask.*median_diameter_band4;
diamHeatmap = diamHeatmap + band5Mask.*median_diameter_band5;
figure('Name', 'diameter heatmap'); set(gcf, 'color', 'w');
imagesc(diamHeatmap); axis image; colormap hot;
caxis([10 18]); xticks([]); yticks([]);
h = colorbar; h.Ticks = [];

% angular alignment
alignHeatmap = zeros(size(band1Mask));
alignHeatmap = alignHeatmap + band1Mask.*median_align_band1;
alignHeatmap = alignHeatmap + band2Mask.*median_align_band2;
alignHeatmap = alignHeatmap + band3Mask.*median_align_band3;
alignHeatmap = alignHeatmap + band4Mask.*median_align_band4;
alignHeatmap = alignHeatmap + band5Mask.*median_align_band5;
% alignHeatmap(alignHeatmap<0) = nan;
figure('Name', 'alignment heatmap'); set(gcf, 'color', 'w');
imagesc(alignHeatmap); axis image; colormap copper;
caxis([0 0.3]); xticks([]); yticks([]);
h = colorbar; h.Ticks = [];

%% get baseline dermis vascular parameters

timePt_list = ["day-3"];
baseline_median_sO2 = zeros(nWound,1);
baseline_median_diameter = zeros(nWound,1);
baseline_median_align = zeros(nWound,1);
baseline_sO2_hetero = zeros(nWound,1);
baseline_median_tortuosity = zeros(nWound,1);
baseline_median_angleChange = zeros(nWound,1);

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)

        structure = load(fullfile(data_base, timePt, 'processed\structure', woundID));
        sO2 = load(fullfile(data_base, timePt, 'processed\sO2', woundID));
        structuralParam = load(fullfile(data_base, timePt, 'processed\paramMaps', woundID));
        tortuosityParam = load(fullfile(data_base, timePt, 'processed\paramMaps\tortuosityMaps', woundID));
        
        structure_map = structure.structure_map_dermis';
        diameter_map = structuralParam.diameterMap;
        align_map = structuralParam.alignMap;
        sO2_map = sO2.sO2_map_dermis';
        structure_map = medfilt2(structure_map,[medfiltSize,medfiltSize]);
        sO2_map = medfilt2(sO2_map,[medfiltSize,medfiltSize]);
        sO2_map(sO2_map>1) = 0;
        sO2_map(sO2_map<0) = 0;
        tortuosity_map = tortuosityParam.tortuosityMap;
        angleChange_map = tortuosityParam.angleChangeMap;
        
        [~,~,sO2Vals_baseline] = find(sO2_map(:));
        median_sO2_baseline = median(sO2Vals_baseline(:), 'omitnan');
        sO2_hetero_baseline = std(sO2Vals_baseline(:), 'omitnan');
        baseline_median_sO2(iWound) = median_sO2_baseline;
        baseline_sO2_hetero(iWound) = sO2_hetero_baseline;
        
        [~,~,diameterVals_baseline] = find(diameter_map(:));
        median_diameter_baseline = median(diameterVals_baseline(:), 'omitnan');
        baseline_median_diameter(iWound) = median_diameter_baseline;
        
        [~,~,alignVals_baseline] = find(align_map(:));
        median_align_baseline = median(alignVals_baseline(:), 'omitnan');
        baseline_median_align(iWound) = median_align_baseline;
        
        [~,~,tortuosityVals_baseline] = find(tortuosity_map(:));
        median_tortuosity_baseline = median(tortuosityVals_baseline(:), 'omitnan');
        baseline_median_tortuosity(iWound) = median_tortuosity_baseline;
        
        [~,~,angleChangeVals_baseline] = find(angleChange_map(:));
        median_angleChange_baseline = median(angleChangeVals_baseline(:), 'omitnan');
        baseline_median_angleChange(iWound) = median_angleChange_baseline;

    end
end

baseline_median_sO2 = baseline_median_sO2 .* 100;
baseline_median_diameter = baseline_median_diameter .* 5;
% baseline_median_align(baseline_median_align<0) = nan;

baseline_sO2_val = median(baseline_median_sO2, 'omitnan');
baseline_sO2_hetero_val = median(baseline_sO2_hetero, 'omitnan');
baseline_diameter_val = median(baseline_median_diameter, 'omitnan');
baseline_align_val = median(baseline_median_align, 'omitnan');
baseline_tortuosity_val = median(baseline_median_tortuosity, 'omitnan');
baseline_angleChange_val = median(baseline_median_angleChange, 'omitnan');

