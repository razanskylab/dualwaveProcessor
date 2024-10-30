%% test script for 3D unmixing
clear, clc; close all;

dwp = dualwaveProcessor();

dwp.data_dir = 'D:\Roy\Experimental data\woundsO2Exp\day-3';
dwp.data_name = 'mouse1_r_5um_dualwave.mat';

dwp.sO2_list = 0 : 0.1 : 1;

dwp.depthStart = 6.5e-3;
dwp.depthEnd = 8.5e-3;

dwp.flagPDResponseScaling = 1;
dwp.flagFiberSplitScaling = 1;
dwp.flagSNRMask = 1;
dwp.flagShowPpeHist = 1;

dwp.fiberSplitRatio = 1.1;
dwp.pdResonsivityRatio = 1.09;
% the product of these two ratios should be around 1.2, which is experimentally calibrated

dwp.amplitudeLowerThresholdPct = 0.04;
dwp.amplitudeUpperThresholdPct = 0.96;

dwp.flagCropFOVBoundary = 1;
dwp.cropXBoundary_mm = 0.5;
dwp.cropYBoundary_mm = 0.5;

dwp.epidermisThickness_mm = 0.012;
dwp.dermisThickness_mm = 0.198;

% take the mode / mean / median of the PPE distributions
dwp.refPpeMethod = 'median';

dwp.Load_Raw_Data();

dwp.Preproc_Basic();

dwp.Correct_Incident_Fluence();

dwp.Load_Tissue_Surface();

dwp.Unmix_3D();

dwp.assignToTissueLayerInfo = 'sO2';
[sO2_map_dermis, sO2_map_hypodermis] = dwp.Assign_To_Tissue_Layer();

dwp.assignToTissueLayerInfo = 'structure';
[structure_map_dermis, structure_map_hypodermis] = dwp.Assign_To_Tissue_Layer();

sO2_map_allDepths = dwp.sO2_map;
structure_map_allDepths = dwp.vesselStructure_map;

% dwp.Render_sO2_3D();
% [structure_map_junction, sO2_map_junction] = dwp.Visualize_Junction_Vasculature();
% dwp.Extract_Tissue_Data();

% visualize structure and sO2 maps with all depths, in dermis and in hypodermis
medfiltKernelSize = 3;
colorbarUB = 60; % [mV]
figure; set(gcf, 'color', 'w');
ax1 = subplot(2,3,1); imagesc(medfilt2(structure_map_allDepths',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax1,bone); colorbar; caxis([0 colorbarUB]); title('Structure image - all depths');
ax2 = subplot(2,3,2); imagesc(medfilt2(structure_map_dermis',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax2,bone); colorbar; caxis([0 colorbarUB]); title('Structure image - dermis');
ax3 = subplot(2,3,3); imagesc(medfilt2(structure_map_hypodermis',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax3,bone); colorbar; caxis([0 colorbarUB]); title('Structure image - hypodermis');
ax4 = subplot(2,3,4); imagesc(medfilt2(sO2_map_allDepths',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax4,turbo); colorbar; caxis([0 1]); title('sO2 image - all depths');
ax5 = subplot(2,3,5); imagesc(medfilt2(sO2_map_dermis',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax5,turbo); colorbar; caxis([0 1]); title('sO2 image - dermis');
ax6 = subplot(2,3,6); imagesc(medfilt2(sO2_map_hypodermis',[medfiltKernelSize,medfiltKernelSize]));
axis image; colormap(ax6,turbo); colorbar; caxis([0 1]); title('sO2 image - hypodermis');

% % save structure and sO2 maps
% mouseInfo = dwp.data_name(1:8);
% save(fullfile(dwp.data_dir, 'processed', 'sO2', mouseInfo), 'sO2_map_allDepths', 'sO2_map_dermis', 'sO2_map_hypodermis');
% save(fullfile(dwp.data_dir, 'processed', 'structure', mouseInfo), 'structure_map_allDepths', 'structure_map_dermis', 'structure_map_hypodermis');

% % show structural and sO2 images side-by-side
% mip532 = squeeze(max(dwp.usVol_532, [], 1));
% mip558 = squeeze(max(dwp.usVol_dye, [], 1));
% sO2_map = dwp.sO2_map;
% sO2_map(sO2_map>1) = 1;
% cbScaleFct = 0.5; % use e.g. 0.5 for better visualization of vasculature
% medfiltKernelSize = 3;
% fUnmix = figure('Name', 'Unmixed image', 'Color', 'w');
% fUnmix.Units = 'normalized';
% fUnmix.OuterPosition = [0 0 1 1];
% cbLimit = dwp.scanSett.sensitivityUs * cbScaleFct;
% ax1 = subplot(131); imagesc(medfilt2(mip532', [medfiltKernelSize,medfiltKernelSize]));
% axis image; colormap(ax1, gray); colorbar; caxis([0 cbLimit]);
% title('532 image used for unmixing');
% ax2 = subplot(132); imagesc(medfilt2(mip558', [medfiltKernelSize,medfiltKernelSize]));
% axis image; colormap(ax2, gray); colorbar; caxis([0 cbLimit]);
% title('Dye image used for unmixing');
% ax3 = subplot(133); imagesc(medfilt2(sO2_map', [medfiltKernelSize,medfiltKernelSize]));
% sO2_LB = min(dwp.sO2_list); sO2_UB = max(dwp.sO2_list);
% axis image; colormap(ax3, turbo); colorbar; caxis([sO2_LB sO2_UB]);
% title('sO2 image');

