%% Script to quantify wound area reduction
clean_all;

data_base = 'D:\Roy\Experimental data\woundsO2Exp';
woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];

pixelSize = 5e-6; % pixel size is 5 um
pixelArea = pixelSize * pixelSize;

woundAreas = zeros(length(woundID_list), length(timePt_list));

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)
        
        bandMasks = load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
        currWoundMask = bandMasks.woundMask;
        
        [counts, binLocs] = imhist(currWoundMask);
        woundArea_mm2 = (counts(1)*pixelArea) * 1e6;

        woundAreas(iWound, iTimePt) = woundArea_mm2;
        
    end
    
end

figure('Name', 'Wound area reduction over time');
set(gcf, 'color', 'w');
H = notBoxPlot(woundAreas, [], 'jitter', 0.3); % convert to [%]
set([H.sdPtch], 'FaceColor', [211,211,211]./255); % [0,0.45,1]
set([H.semPtch], 'FaceColor', [64,64,64]./255); % [1,0.55,0]
set([H.mu],'color','w');
set([H.data],'LineWidth',1);
ylim([0 20]); grid on;
xticklabels({'day 4', 'day 7', 'day 10'});
tickLabels = get(gca,'XTickLabel');
set(gca, 'XTickLabel', tickLabels, 'fontsize', 14);
ax = gca; ax.LineWidth = 2;

[h,p] = ttest(woundAreas(:,1), woundAreas(:,2))
[h,p] = ttest(woundAreas(:,1), woundAreas(:,3))
[h,p] = ttest(woundAreas(:,2), woundAreas(:,3))


