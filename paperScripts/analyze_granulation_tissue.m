%% Script to quantify the vessel parameters in the granulation tissue area
clean_all;

data_base = 'C:\data_temp\woundExp';
woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];

pixelSize = 5e-6; % pixel size is 5 um
pixelArea = pixelSize * pixelSize;
woundDiam_day0 = 5e-3;
woundRadiIdx_day0 = round(woundDiam_day0/2/pixelSize);
% make nTheta points along the circular boundary of the wound at day0
nTheta = 1000;
deltaTheta = (2*pi) / nTheta;
theta = 0 : deltaTheta : 2*pi;

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)

        % load structural image
        load(fullfile(data_base, timePt, 'processed/structure', woundID));
        nPixels = size(structure_map_allDepths,1);
        figure; imagesc(structure_map_allDepths'); axis image;
        colormap bone; colorbar; caxis([0 80]); hold on;
        title(strcat(woundID, '_', timePt));

        % get band masks
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        woundMask = bandMasks.woundMask;
        [counts, binLocs] = imhist(woundMask);
        woundArea_mm2 = (counts(1)*pixelArea) * 1e6;

        % get wound boundary
        currBoundary = ReadImageJROI(fullfile(data_base, timePt, 'processed\forMasking', strcat(woundID, '.roi')));
        boundaryCoordinates = currBoundary.mnCoordinates;
        geoCenter = round(mean(boundaryCoordinates, 1));
        scatter(boundaryCoordinates(:,1),boundaryCoordinates(:,2),'r'); hold on;

        % get the circular wound boundary at day0 with the center at the geocenter
        % of current wound boundary
        xCoorCircle = cos(theta) .* woundRadiIdx_day0 + geoCenter(1);
        yCoorCircle = sin(theta) .* woundRadiIdx_day0 + geoCenter(2);
        % take care of boundary condition
        xCoorCircle(xCoorCircle<1) = 1;
        xCoorCircle(xCoorCircle>nPixels) = nPixels;
        yCoorCircle(yCoorCircle<1) = 1;
        yCoorCircle(yCoorCircle>nPixels) = nPixels;
        scatter(xCoorCircle, yCoorCircle, 'g'); hold off;
        drawnow;

        % get circulate wound boundary at day0 and save
        [xq, yq] = meshgrid(1:nPixels, 1:nPixels);
        xq = xq(:); yq = yq(:);
        woundMask_day0 = inpolygon(xq, yq, xCoorCircle, yCoorCircle);
        woundMask_day0 = reshape(~woundMask_day0,[nPixels,nPixels]);
        
        save(fullfile(data_base, timePt, 'processed\woundMaskDay0', woundID), 'woundMask_day0');
        
    end
    
end

%% quantify vessel parameters within granulation tissue areas

nWound = length(woundID_list);
nTimePt = length(timePt_list);

sO2_stat = zeros(nWound, nTimePt);
sO2_hetero_stat = zeros(nWound, nTimePt);
diameter_stat = zeros(nWound, nTimePt);
diameter_hetero_stat = zeros(nWound, nTimePt);
align_stat = zeros(nWound, nTimePt);
tortuosity_stat = zeros(nWound, nTimePt);
angleChange_stat = zeros(nWound, nTimePt);

medfiltSize = 3;

for iWound = 1 : nWound
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : nTimePt
        timePt = timePt_list(iTimePt)
        
        % get corresponding hand-drawn wound mask
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        woundMask = bandMasks.woundMask;
        
        % get corresponding original circular mask
        originalMask = load(fullfile(data_base, timePt, 'processed\woundMaskDay0', woundID));
        woundMask_day0 = originalMask.woundMask_day0;
        
        % get the granulation tissue area
        granulationMask = (~woundMask_day0) .* woundMask;
%         save(fullfile(data_base, timePt, 'processed\granulationMask', woundID), 'granulationMask');
        
        % extract vessel parameters within the granulation area
        structure = load(fullfile(data_base, timePt, 'processed\structure', woundID));
        sO2 = load(fullfile(data_base, timePt, 'processed\sO2', woundID));
        structuralParam = load(fullfile(data_base, timePt, 'processed\paramMaps', woundID));
        newAlignment = load(fullfile(data_base, timePt, 'processed\paramMaps\newAlignMaps', woundID));
        tortuosityParam = load(fullfile(data_base, timePt, 'processed\paramMaps\tortuosityMaps', woundID));
        
        % structural information
        structure_map = structure.structure_map_dermis';
        structure_map = structure_map .* granulationMask;
        structure_map = medfilt2(structure_map, [medfiltSize,medfiltSize]);

        % oxygenation information
        sO2_map = sO2.sO2_map_dermis';
        sO2_map = sO2_map .* granulationMask;
        sO2_map = medfilt2(sO2_map, [medfiltSize,medfiltSize]);
        % remove nonsense sO2 values
        sO2_map(sO2_map>1) = 0;
        sO2_map(sO2_map<0) = 0;
        
        % diameter information
        diameter_map = structuralParam.diameterMap;
        diameter_map = diameter_map .* granulationMask;
        
        % angular alignment information
        align_map = newAlignment.alignMap;
        align_map = align_map .* granulationMask;
        
        % tortuosity information
        tortuosity_map = tortuosityParam.tortuosityMap;
        tortuosity_map = tortuosity_map .* granulationMask;
        
        % angular change information
        angleChange_map = tortuosityParam.angleChangeMap;
        angleChange_map = angleChange_map .* granulationMask;
        
        [~,~,sO2Vals] = find(sO2_map(:));
        median_sO2 = median(sO2Vals(:), 'omitnan');
        std_sO2 = std(sO2Vals(:), 'omitnan');
        [~,~,diameterVals] = find(diameter_map(:));
        median_diameter = median(diameterVals(:), 'omitnan');
        std_diameter = std(diameterVals(:), 'omitnan');
        [~,~,alignVals] = find(align_map(:));
        median_align = median(alignVals(:), 'omitnan');
        [~,~,tortuosityVals] = find(tortuosity_map(:));
        median_tortuosity = median(tortuosityVals(:), 'omitnan');
        [~,~,angleChangeVals] = find(angleChange_map(:));
        median_angleChange = median(angleChangeVals(:), 'omitnan');
        
        sO2_stat(iWound, iTimePt) = median_sO2;
        sO2_hetero_stat(iWound, iTimePt) = std_sO2;
        diameter_stat(iWound, iTimePt) = median_diameter;
        diameter_hetero_stat(iWound, iTimePt) = std_diameter;
        align_stat(iWound, iTimePt) = median_align;
        tortuosity_stat(iWound, iTimePt) = median_tortuosity;
        angleChange_stat(iWound, iTimePt) = median_angleChange;
        
    end
    
end

sO2_stat = sO2_stat .* 100; % [%]
sO2_hetero_stat = sO2_hetero_stat .* 100; % [%]
diameter_stat = diameter_stat .* 5; % [um]
diameter_hetero_stat = diameter_hetero_stat .* 5; % [um]

save('C:\data_temp\woundExp\granuStat.mat', 'sO2_stat', 'diameter_stat', ...
    'align_stat', 'tortuosity_stat', 'angleChange_stat');

%% generate box plot and quantify significance

figure;
subplot(131); boxplot(sO2_stat); title('sO2 changes in granulation area');
subplot(132); boxplot(diameter_stat); title('Diameter changes in granulation area');
subplot(133); boxplot(align_stat); title('Alignment changes in granulation area');

figure('Name', 'sO2 temporal changes within granulation tissue');
set(gcf, 'color', 'w');
H = notBoxPlot(sO2_stat, [], 'jitter', 0.3);
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([30 80]); grid on; hold on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_sO2_val, '--', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

figure('Name', 'Diameter temporal changes within granulation tissue');
set(gcf, 'color', 'w');
H = notBoxPlot(diameter_stat, [], 'jitter', 0.3); % convert to [%]
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([10 18]); grid on; hold on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_diam_val, '--', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

figure('Name', 'Angular alignment temporal changes within granulation tissue');
set(gcf, 'color', 'w');
H = notBoxPlot(align_stat, [], 'jitter', 0.3); % convert to [%]
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([0 0.5]); grid on;
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

figure('Name', 'tortuosity temporal changes within granulation tissue');
set(gcf, 'color', 'w');
H = notBoxPlot(tortuosity_stat, [], 'jitter', 0.3); % convert to [%]
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([1.025 1.05]); grid on;
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

figure('Name', 'angleChange temporal changes within granulation tissue');
set(gcf, 'color', 'w');
H = notBoxPlot(angleChange_stat, [], 'jitter', 0.3); % convert to [%]
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([6 11]); grid on;
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

% use paired ttest for longitudinal data
[h,p] = ttest(sO2_stat(:,1), sO2_stat(:,2))
[h,p] = ttest(sO2_stat(:,1), sO2_stat(:,3))
[h,p] = ttest(sO2_stat(:,2), sO2_stat(:,3))

[h,p] = ttest(diameter_stat(:,1), diameter_stat(:,2))
[h,p] = ttest(diameter_stat(:,1), diameter_stat(:,3))
[h,p] = ttest(diameter_stat(:,2), diameter_stat(:,3))

[h,p] = ttest(align_stat(:,1), align_stat(:,2))
[h,p] = ttest(align_stat(:,1), align_stat(:,3))
[h,p] = ttest(align_stat(:,2), align_stat(:,3))

