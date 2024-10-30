%% script to generate overlapping histogram

data_base = 'D:\Experimental data\Skull data';
strain = 'Nude';
%strain = 'Black6';
%strain = 'CD1';

fileName = strcat('VA_', strain, 'W4toW16_Frangi_040221.mat');

mouse_list = ["M01", "M02", "M03", "M10", "M20"];

% load stat.
stat = load(fullfile(data_base, strain, fileName));

% set colors
c1 = [255, 0, 90] ./ 255;
c2 = [255, 165, 0] ./ 255;
c3 = [0, 255, 165] ./ 255;
c4 = [0, 90, 255] ./ 255;

% set bin number and limit
x_start = 1; x_end = 15; 
x_step = 0.25;
xLimit = x_start : x_step : x_end;
normalization = 'count';
doFit = 0;
nBins = length(xLimit);

for i_mouse = 1 : length(mouse_list)

    % plot overlapping histogram - Nude
    figure(1); set(gcf, 'color', 'w');
    hw4 = pretty_hist(stat.(strcat(strain,'W04',mouse_list(i_mouse))).vesDiameter, c1, nBins, normalization, doFit, xLimit);
    hold on;
    hw8 = pretty_hist(stat.(strcat(strain,'W08',mouse_list(i_mouse))).vesDiameter, c2, nBins, normalization, doFit, xLimit);
    hold on;
    hw12 = pretty_hist(stat.(strcat(strain,'W12',mouse_list(i_mouse))).vesDiameter, c3, nBins, normalization, doFit, xLimit);
    hold on;
    hw16 = pretty_hist(stat.(strcat(strain,'W16',mouse_list(i_mouse))).vesDiameter, c4, nBins, normalization, doFit, xLimit);
    hold off; grid off;
    legend([hw4,hw8,hw12,hw16],{'Week 4','Week 8','Week 12','Week 16'}, 'FontSize', 20);
    xlim([x_start, x_end]);
    ylim([0 250])
    % xticks([]); 
    yticks([]);
    export_name = strcat(strain, '_', mouse_list(i_mouse));
    export_fig(export_name, '-jpg', '-r300');
    
end

% % plot overlapping histogram - Black6
% figure; set(gcf, 'color', 'w');
% hw4 = pretty_hist(stat.Black6W04M02.vesDiameter, c1, nBins, normalization, doFit, xLimit);
% hold on;
% hw8 = pretty_hist(stat.Black6W08M02.vesDiameter, c2, nBins, normalization, doFit, xLimit);
% hold on;
% hw12 = pretty_hist(stat.Black6W12M02.vesDiameter, c3, nBins, normalization, doFit, xLimit);
% hold on;
% hw16 = pretty_hist(stat.Black6W16M02.vesDiameter, c4, nBins, normalization, doFit, xLimit);
% hold off; grid off;
% legend([hw4,hw8,hw12,hw16],{'Week 4','Week 8','Week 12','Week 16'}, 'FontSize', 14);
% xlim([x_start, x_end]);
% ylim([0 350])
% % xticks([]); 
% yticks([]);

% % plot overlapping histogram - CD1
% figure; set(gcf, 'color', 'w');
% hw4 = pretty_hist(stat.CD1W04M02.vesDiameter, c1, nBins, normalization, doFit, xLimit);
% hold on;
% hw8 = pretty_hist(stat.CD1W08M02.vesDiameter, c2, nBins, normalization, doFit, xLimit);
% hold on;
% hw12 = pretty_hist(stat.CD1W12M02.vesDiameter, c3, nBins, normalization, doFit, xLimit);
% hold on;
% hw16 = pretty_hist(stat.CD1W16M02.vesDiameter, c4, nBins, normalization, doFit, xLimit);
% hold off; grid off;
% legend([hw4,hw8,hw12,hw16],{'Week 4','Week 8','Week 12','Week 16'}, 'FontSize', 14);
% xlim([x_start, x_end]);
% ylim([0 350])
% % xticks([]); 
% yticks([]);

