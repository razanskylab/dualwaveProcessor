%% generate new set of box plots with the new boxchart function
% for wound area reduction
clear, clc; close all;

load('woundAreaStat.mat');

pointJitter = 0.01; % lateral jittering of sample points for better visibility
scatterLineWidth = 0.76;
scatterMarkerFaceColor = [0,0.7,0.6];
scaterMarkerEdgeColor = [0,0,0];
boxWidth = 0.4;
boxFaceColor = [0.1,0.1,0.1];
boxLineWidth = 1;

[nWound, nTimePt] = size(woundAreas);

oriWoundDiam = 5; % [mm]
oriWoundRadius = oriWoundDiam / 2;
oriWoundSize = pi*(oriWoundRadius^2);

figure('Name', 'Wound area reduction over time'); set(gcf, 'color', 'w');
H = boxchart(woundAreas); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*oriWoundSize, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
H.MarkerStyle = 'none'; % do NOT show the outlier as an additional circle
yticks([0 : 5 : 20]);
ax = gca; ax.LineWidth = 2;

%% generate new set of box plots with the new boxchart function
% For all spatial changes
clear, clc; close all;

load('summaryStat.mat');
[nWound, nTimePt, nBand] = size(sO2_stat);

pointJitter = 0.025; % lateral jittering of sample points for better visibility
scatterLineWidth = 0.76;
scatterMarkerFaceColor = [0,0.7,0.6];
scaterMarkerEdgeColor = [0,0,0];
boxWidth = 0.3;
boxFaceColor = [0.1,0.1,0.1];
boxLineWidth = 1;

pixelSize = 5; % [um]

% sO2 spatial changes at day 4
sO2_spatial_day4 = squeeze(sO2_stat(:,1,:));
baseline_sO2 = mean(baseline_median_sO2);
figure('Name', 'sO2 Spatial Change at day 4'); set(gcf, 'color', 'w');
H = boxchart(sO2_spatial_day4); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nBand+1, ones(1,nBand+2).*baseline_sO2, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
% title('sO2 spatial changes at day 4');
ylim([30 80]); yticks([30 : 10 : 80]);

% diameter spatial changes at day 4
diameter_spatial_day4 = squeeze(diameter_stat(:,1,:));
baseline_diameter = mean(baseline_median_diameter);
figure('Name', 'Diameter Spatial Change at day 4'); set(gcf, 'color', 'w');
H = boxchart(diameter_spatial_day4); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nBand+1, ones(1,nBand+2).*baseline_diameter, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
% title('diameter spatial changes at day 4');
ylim([10 18]); yticks([10 : 2 : 18]);

% tortuosity spatial changes at day 4
tortuosity_spatial_day4 = squeeze(angleChange_stat(:,1,:) ./ pixelSize);
baseline_tortuosity = mean(baseline_median_angleChange./pixelSize);
figure('Name', 'Tortuosity Spatial Change at day 4'); set(gcf, 'color', 'w');
H = boxchart(tortuosity_spatial_day4); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nBand+1, ones(1,nBand+2).*baseline_tortuosity, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
% title('tortuosity spatial changes at day 4');
ylim([1 2.5]); yticks([1 : 0.5 : 2.5]);

% % angular alignment spatial changes at day 4
% align_spatial_day4 = squeeze(align_stat(:,1,:));
% baseline_align = mean(baseline_median_align);
% figure('Name', 'Tortuosity Spatial Change at day 4'); set(gcf, 'color', 'w');
% H = boxchart(align_spatial_day4); hold on;
% Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
% scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
%     'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
% hold on; grid on;
% plot(0:nBand+1, ones(1,nBand+2).*baseline_align, '--', 'LineWidth', 2, ...
%     'Color', [0.8500 0.3250 0.0980]); hold off;
% H.BoxWidth = boxWidth;
% H.BoxFaceColor = boxFaceColor;
% H.LineWidth = boxLineWidth;
% ax = gca; ax.LineWidth = 2;
% title('alignment spatial changes at day 4');

%% For temporal changes in the healing front (i.e. band 1)
clear, clc; close all;

load('summaryStat.mat');
[nWound, nTimePt, nBand] = size(sO2_stat);

pointJitter = 0.03; % lateral jittering of sample points for better visibility
scatterLineWidth = 0.76;
scatterMarkerFaceColor = [0,0.7,0.6];
scaterMarkerEdgeColor = [0,0,0];
boxWidth = 0.3;
boxFaceColor = [0.1,0.1,0.1];
boxLineWidth = 1;

pixelSize = 5;

% sO2 temporal changes in band 1
sO2_temporal_band1 = squeeze(sO2_stat(:,:,1)); % [nWound, nTimePt]
baseline_sO2 = mean(baseline_median_sO2);
figure('Name', 'sO2 Temporal Change in band 1'); set(gcf, 'color', 'w');
H = boxchart(sO2_temporal_band1); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_sO2, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([30 80]); yticks([30 : 10 : 80]);

% diameter temporal changes in band 1
diameter_temporal_band1 = squeeze(diameter_stat(:,:,1)); % [nWound, nTimePt]
baseline_diameter = mean(baseline_median_diameter);
figure('Name', 'Diameter Temporal Change in band 1'); set(gcf, 'color', 'w');
H = boxchart(diameter_temporal_band1); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_diameter, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([10 18]); yticks([10 : 2 : 18]);

% turtuosity temporal changes in band 1
tortuosity_temporal_band1 = squeeze(angleChange_stat(:,:,1)./pixelSize); % [nWound, nTimePt]
baseline_tortuosity = mean(baseline_median_angleChange./pixelSize);
figure('Name', 'Tortuosity Temporal Change in band 1'); set(gcf, 'color', 'w');
H = boxchart(tortuosity_temporal_band1); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_tortuosity, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([1 2.5]); yticks([1 : 0.5 : 2.5]);

% % angular alignment temporal changes in band 1
% align_temporal_band1 = squeeze(align_stat(:,:,1)); % [nWound, nTimePt]
% figure('Name', 'Angular alignment Temporal Change in band 1'); set(gcf, 'color', 'w');
% H = boxchart(align_temporal_band1); hold on;
% Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
% scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
%     'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
% hold off;
% H.BoxWidth = boxWidth;
% H.BoxFaceColor = boxFaceColor;
% H.LineWidth = boxLineWidth;
% ax = gca; ax.LineWidth = 2;

%% For temporal changes in the granulation tissue
clear, clc; close all;

load('granuStat.mat');
[nWound, nTimePt] = size(sO2_stat);

pointJitter = 0.03; % lateral jittering of sample points for better visibility
scatterLineWidth = 0.76;
scatterMarkerFaceColor = [0,0.7,0.6];
scaterMarkerEdgeColor = [0,0,0];
boxWidth = 0.3;
boxFaceColor = [0.1,0.1,0.1];
boxLineWidth = 1;

pixelSize = 5;

bandStat = matfile('summaryStat.mat');
baseline_sO2 = mean(bandStat.baseline_median_sO2);
baseline_diameter = mean(bandStat.baseline_median_diameter);
baseline_tortuosity = mean(bandStat.baseline_median_angleChange./pixelSize);

% sO2 temporal changes in granulation tissue
figure('Name', 'sO2 Temporal Change in granulation tissue'); set(gcf, 'color', 'w');
H = boxchart(sO2_stat); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_sO2, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([30 80]); yticks([30 : 10 : 80]);

% diameter temporal changes in granulation tissue
figure('Name', 'diameter Temporal Change in granulation tissue'); set(gcf, 'color', 'w');
H = boxchart(diameter_stat); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_diameter, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([10 18]); yticks([10 : 2 : 18]);
H.MarkerStyle = 'none';

% tortuosity temporal changes in granulation tissue
figure('Name', 'tortuosity Temporal Change in granulation tissue'); set(gcf, 'color', 'w');
H = boxchart(angleChange_stat./pixelSize); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_tortuosity, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([1 2.5]); yticks([1 : 0.5 : 2.5]);

%% For new supplementary figures on diameter heterogeneity and angular alignment
clear, clc; close all;

load('supplementaryDiamHeteroAlign_stat.mat');
diameterHetero_stat = diameterHetero_stat .* 5;

[nWound, nTimePt, nBand] = size(diameterHetero_stat);

pointJitter = 0.01; % lateral jittering of sample points for better visibility
scatterLineWidth = 0.76;
scatterMarkerFaceColor = [0,0.7,0.6];
scaterMarkerEdgeColor = [0,0,0];
boxWidth = 0.4;
boxFaceColor = [0.1,0.1,0.1];
boxLineWidth = 1;

% diameter spatial changes at day 4
diameterHetero_spatial_day4 = squeeze(diameterHetero_stat(:,1,:));
baseline_diameterHetero = mean(baseline_diameter_hetero);
figure('Name', 'Diameter Heterogeneity Spatial Change at day 4'); set(gcf, 'color', 'w');
H = boxchart(diameterHetero_spatial_day4); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nBand+1, ones(1,nBand+2).*baseline_diameterHetero, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
% title('diameter spatial changes at day 4');
ylim([2.5, 17.5]); yticks([2.5 : 5 : 17.5]);

% diameter temporal changes in band 1
diameterHetero_temporal_band1 = squeeze(diameterHetero_stat(:,:,1)); % [nWound, nTimePt]
figure('Name', 'Diameter Heterogeneity Temporal Change in band 1'); set(gcf, 'color', 'w');
H = boxchart(diameterHetero_temporal_band1); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold on; grid on;
plot(0:nTimePt+1, ones(1,nTimePt+2).*baseline_diameterHetero, '--', 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980]); hold off;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([2.5, 17.5]); yticks([2.5 : 5 : 17.5]);

% angular alignment spatial changes at day 10
align_spatial_day10 = squeeze(align_stat(:,3,:));
figure('Name', 'Angular alignment Spatial Change at day 10'); set(gcf, 'color', 'w');
H = boxchart(align_spatial_day10); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nBand).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold off; grid on;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
% title('diameter spatial changes at day 4');
ylim([0, 0.5]); yticks([0 : 0.1 : 0.5]);

% angular alignment temporal changes in band 3
align_temporal_band3 = squeeze(align_stat(:,:,3)); % [nWound, nTimePt]
figure('Name', 'Angular alignment Temporal Change in band 3'); set(gcf, 'color', 'w');
H = boxchart(align_temporal_band3); hold on;
Xjittered = repmat(double(H.XData),nWound,1) + randn(nWound,nTimePt).*pointJitter;
scatter(Xjittered, H.YData, 'filled', 'LineWidth', scatterLineWidth, ...
    'MarkerFaceColor', scatterMarkerFaceColor, 'MarkerEdgeColor',scaterMarkerEdgeColor);
hold off; grid on;
H.BoxWidth = boxWidth;
H.BoxFaceColor = boxFaceColor;
H.LineWidth = boxLineWidth;
ax = gca; ax.LineWidth = 2;
ylim([0, 0.5]); yticks([0 : 0.1 : 0.5]);

%% plot updated line graphs for supplementary figures

clear, clc; close all;
load('supplementaryDiamHeteroAlign_stat.mat');
[nSample, nTimePt, nBand] = size(diameterHetero_stat);
pixelSize = 5;
diameterHetero_stat = diameterHetero_stat .* pixelSize;
scatterSize = 20;

figure('Name', 'Summary line plot for diameter heterogeneities changes'); set(gcf, 'color', 'w');
summary_diameterHetero = squeeze(median(diameterHetero_stat,1));
curve_day4 = summary_diameterHetero(1,:);
curve_day7 = summary_diameterHetero(2,:);
curve_day10 = summary_diameterHetero(3,:);
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(diameterHetero_stat(:,1,:));
YData = YData(:);
% plot(curve_day4, 'LineWidth', 3, 'Color', [0 0.4470 0.7410]); hold on;
moe = squeeze(diameterHetero_stat(:,1,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day4, moe, 'LineWidth', 3, 'Color', [0 0.4470 0.7410]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'none'); hold on;
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(diameterHetero_stat(:,2,:));
YData = YData(:);
% plot(curve_day7, 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250]); hold on;
moe = squeeze(diameterHetero_stat(:,2,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day7, moe, 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', 'none'); hold on;
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(diameterHetero_stat(:,3,:));
YData = YData(:);
% plot(curve_day10, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880]); hold on;
moe = squeeze(diameterHetero_stat(:,3,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day10, moe, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', 'none'); hold on;
plot(ones(size(curve_day4)).*baseline_diameter_hetero_val, '--', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]); hold off;
ylim([2.5, 17.5]); xticks([1:nBand]);
% legend('day 4', 'day 7', 'day 10', 'baseline', 'FontSize', 14);
xticks([]); yticks([]); set(gca, 'Visible', 'Off');

figure('Name', 'Summary line plot for angular alignment changes'); set(gcf, 'color', 'w');
summary_align = squeeze(median(align_stat,1));
curve_day4 = summary_align(1,:);
curve_day7 = summary_align(2,:);
curve_day10 = summary_align(3,:);
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(align_stat(:,1,:));
YData = YData(:);
% plot(curve_day4, 'LineWidth', 3, 'Color', [0 0.4470 0.7410]); hold on;
moe = squeeze(align_stat(:,1,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day4, moe, 'LineWidth', 3, 'Color', [0 0.4470 0.7410]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'none'); hold on;
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(align_stat(:,2,:));
YData = YData(:);
% plot(curve_day7, 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250]); hold on;
moe = squeeze(align_stat(:,2,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day7, moe, 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', 'none'); hold on;
XData = ones(nSample,nBand) .* (1:nBand);
XData = XData(:);
YData = squeeze(align_stat(:,3,:));
YData = YData(:);
% plot(curve_day10, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880]); hold on;
moe = squeeze(align_stat(:,3,:)); % margin of error
moe = 1.96 .* (std(moe,1)./sqrt(nSample)); % corresponds to 95% confidence interval
errorbar(curve_day10, moe, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880]); hold on;
scatter(XData, YData, scatterSize, 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', 'none'); hold off;
ylim([0, 0.5]); xticks([1:nBand]);
% legend('day 4', 'day 7', 'day 10', 'baseline', 'FontSize', 14);
xticks([]); yticks([]); set(gca, 'Visible', 'Off');
