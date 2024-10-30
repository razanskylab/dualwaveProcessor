%% Script to segment multi-layered vessels in OA dataset based on estimated skin surface
clean_all;

data_base = 'D:\Roy\Experimental data\woundsO2Exp\day-3';
data_name = 'mouse1_l_5um_dualwave.mat';

dwp = dualwaveProcessor();
dwp.data_dir = data_base;
dwp.data_name = data_name;

dwp.depthStart = 6e-3;
dwp.depthEnd = 8.5e-3;

dwp.flagCropFOVBoundary = 1;
dwp.cropXBoundary_mm = 0.5;
dwp.cropYBoundary_mm = 0.5;

dwp.Load_Raw_Data();
dwp.Preproc_Basic();

%% Load estimated skin surface and show overlay on some b-scans

% NOTE: watch out the value here!
% depthCorrection_mm = -0.018; % this is the assumed epidermis thickness in 3D unmixing
depthCorrection_mm = 0;

mouseInfo = data_name(1:8);
load(fullfile(data_base, 'skinCropDepths', mouseInfo));
[nAlines_us, nBscan_us] = size(skinCropDepths_mm);
[~, nAlines_oa, nBscan_oa] = size(dwp.usVol_532);

skinCropDepths_mm = skinCropDepths_mm + depthCorrection_mm;

% interpolate US scanning grid to match finer OA scanning grid
stepSize_mm_us = 10e-3;
stepSize_mm_oa = 5e-3;
xVec_mm_us = (0:nAlines_us-1) .* stepSize_mm_us;
yVec_mm_us = (0:nBscan_us-1) .* stepSize_mm_us;
xVec_mm_oa = (0:nAlines_oa-1) .* stepSize_mm_oa;
yVec_mm_oa = (0:nBscan_oa-1) .* stepSize_mm_oa;
usGrid{1} = xVec_mm_us;
usGrid{2} = yVec_mm_us;
oaGrid{1} = xVec_mm_oa;
oaGrid{2} = yVec_mm_oa;
interpF = griddedInterpolant(usGrid, skinCropDepths_mm, 'linear', 'linear');
skinCropDepths_mm_interp = interpF(oaGrid);

bscanInterval_mm = 0.2;
dermisThickness_mm = 0.204;
fullThickness_mm = 0.5;
bscanInterval_idx = bscanInterval_mm / stepSize_mm_oa;
zVec_mm = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
zVec_mm = zVec_mm(dwp.zCropIdx_us(1,1):dwp.zCropIdx_us(1,2));
for iBscan_oa = 1 : nBscan_oa
    currBscan = dwp.usVol_532(:, :, iBscan_oa);

    % show some B-scans together with the cropping profile
    if ~mod(iBscan_oa,bscanInterval_idx)
        fBscan = figure('Name', 'B-scan', 'Color', 'w');
        fBscan.Units = 'normalized';
        fBscan.OuterPosition = [0 0 1 1];
        imagesc(xVec_mm_oa, zVec_mm, currBscan);
        colormap redblue; colorbar; axis tight equal;
        caxis([-40 40]);
        xlabel('[mm]'); ylabel('[mm]');
        title(sprintf('%i-th B-scan', iBscan_oa));
        hold on;
        currSkinSurface_mm = skinCropDepths_mm_interp(:,iBscan_oa);
        plot(xVec_mm_oa, currSkinSurface_mm, 'g', 'LineWidth', 2);

        % also show different layers
        hold on;
        plot(xVec_mm_oa, currSkinSurface_mm+dermisThickness_mm, 'g--', 'LineWidth', 2);
        
        hold on;
        plot(xVec_mm_oa, currSkinSurface_mm+fullThickness_mm, 'k--', 'LineWidth', 2);
        
        hold off;
        legend('Skin surface', 'dermis-hypodermis junction', 'End of skin');
        
        waitforbuttonpress;
        close;
    end
    
end
% update saved skin surface profile with double checked profile
save(fullfile(data_base, 'skinCropDepths', mouseInfo), 'skinCropDepths_mm');

%% Visualize multi-layered vessels starting from the estimated skin surface

epidermisDepth_mm = 0.006;
dermisDepth_mm = 0.198;
hypodermisDepth_mm = 0.198;

% usVol has size: [nt, nAlines, nBscan], crop depths has size: [nAlines, nBscan]
skinCropIndices_interp = zeros(size(skinCropDepths_mm_interp));
for iBscan = 1 : nBscan_oa
    for iAline = 1 : nAlines_oa
        
        [~,currCropIdx] = min(abs(zVec_mm - skinCropDepths_mm_interp(iAline, iBscan)));
        skinCropIndices_interp(iAline, iBscan) = currCropIdx;
        
    end
end

samplingDist_mm = (1/dwp.scanSett.samplingFreq)*dwp.scanSett.SOS*1e3;
nPixelEpidermis = round(epidermisDepth_mm / samplingDist_mm);
if mod(nPixelEpidermis,2)
    % make sure this is an even number
    nPixelEpidermis = nPixelEpidermis + 1;
end
nPixelDermis = ceil(dermisDepth_mm / samplingDist_mm);
nPixelHypodermis = ceil(hypodermisDepth_mm / samplingDist_mm);

epidermisVol = zeros(nPixelEpidermis, nAlines_oa, nBscan_oa);
dermisVol = zeros(nPixelDermis, nAlines_oa, nBscan_oa);
hypodermisVol = zeros(nPixelHypodermis, nAlines_oa, nBscan_oa);

nz = length(zVec_mm);

for iBscan = 1 : nBscan_oa
    for iAline = 1 : nAlines_oa
        
        currALine = dwp.usVol_532(:, iAline, iBscan);
        currSkinStartIdx = skinCropIndices_interp(iAline, iBscan);
        epidermisStartIdx = currSkinStartIdx - (nPixelEpidermis/2) + 1;
        epidermisEndIdx = currSkinStartIdx + (nPixelEpidermis/2);
        dermisStartIdx = epidermisEndIdx + 1;
        dermisEndIdx = dermisStartIdx + nPixelDermis - 1;
        hypodermisStartIdx = dermisEndIdx + 1;
        hypodermisEndIdx = hypodermisStartIdx + nPixelHypodermis - 1;
        
        endIndices = [epidermisEndIdx,dermisEndIdx,hypodermisEndIdx];
        if any(endIndices>nz)
            continue;
        else
            epidermisVol(:,iAline,iBscan) = currALine(epidermisStartIdx:epidermisEndIdx);
            dermisVol(:,iAline,iBscan) = currALine(dermisStartIdx:dermisEndIdx);
            hypodermisVol(:,iAline,iBscan) = currALine(hypodermisStartIdx:hypodermisEndIdx);
        end
        
    end
end

colorUB = 100;
mipAllDepth = squeeze(max(dwp.usVol_532,[],1))';
mipEpidermis = squeeze(max(epidermisVol,[],1))';
mipDermis = squeeze(max(dermisVol,[],1))';
mipHypodermis = squeeze(max(hypodermisVol,[],1))';

fH(1) = figure('Name', 'All depths'); set(gcf, 'color', 'w');
imagesc(mipAllDepth); axis image; colormap inferno; h = colorbar; h.Ticks = [];
xticks([]); yticks([]); caxis([0 colorUB]); title('MIP over all depths');
fH(2) = figure('Name', 'Epidermis'); set(gcf, 'color', 'w');
imagesc(mipEpidermis); axis image; colormap inferno; h = colorbar; h.Ticks = [];
xticks([]); yticks([]); caxis([0 colorUB]); title('MIP over epidermis');
fH(3) = figure('Name', 'Dermis'); set(gcf, 'color', 'w');
imagesc(mipDermis); axis image; colormap inferno; h = colorbar; h.Ticks = [];
xticks([]); yticks([]); caxis([0 colorUB]); title('MIP over dermis');
fH(4) = figure('Name', 'Hypodermis'); set(gcf, 'color', 'w');
imagesc(mipHypodermis); axis image; colormap inferno; h = colorbar; h.Ticks = [];
xticks([]); yticks([]); caxis([0 colorUB]); title('MIP over hypodermis');

% save segmented maps
saveDir = fullfile(data_base, 'processed\structure');
save(fullfile(saveDir, mouseInfo), 'mipAllDepth', 'mipEpidermis', 'mipDermis', 'mipHypodermis');

%% show sO2 maps over time
data_base = 'C:\data_temp\woundExp';
medfiltSize = 3;

woundID = 'mouse4_r';

timePt = 'day-3';
load(fullfile(data_base, timePt, 'processed', 'sO2', woundID));
figure('Name', woundID); set(gcf, 'color', 'w');
ax1 = subplot(3,4,1); imagesc(medfilt2(sO2_map_allDepths',[medfiltSize,medfiltSize]));
axis image; colormap(ax1, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax5 = subplot(3,4,5); imagesc(medfilt2(sO2_map_dermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax5, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax9 = subplot(3,4,9); imagesc(medfilt2(sO2_map_hypodermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax9, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);

timePt = 'day4';
load(fullfile(data_base, timePt, 'processed', 'sO2', woundID));
ax2 = subplot(3,4,2); imagesc(medfilt2(sO2_map_allDepths',[medfiltSize,medfiltSize]));
axis image; colormap(ax2, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax6 = subplot(3,4,6); imagesc(medfilt2(sO2_map_dermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax6, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax10 = subplot(3,4,10); imagesc(medfilt2(sO2_map_hypodermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax10, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);

timePt = 'day7';
load(fullfile(data_base, timePt, 'processed', 'sO2', woundID));
ax3 = subplot(3,4,3); imagesc(medfilt2(sO2_map_allDepths',[medfiltSize,medfiltSize]));
axis image; colormap(ax3, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax7 = subplot(3,4,7); imagesc(medfilt2(sO2_map_dermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax7, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);
ax11 = subplot(3,4,11); imagesc(medfilt2(sO2_map_hypodermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax11, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt);

timePt = 'day10';
load(fullfile(data_base, timePt, 'processed', 'sO2', woundID));
ax4 = subplot(3,4,4); imagesc(medfilt2(sO2_map_allDepths',[medfiltSize,medfiltSize]));
axis image; colormap(ax4, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt); colorbar;
ax8 = subplot(3,4,8); imagesc(medfilt2(sO2_map_dermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax8, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt); colorbar;
ax12 = subplot(3,4,12); imagesc(medfilt2(sO2_map_hypodermis',[medfiltSize,medfiltSize]));
axis image; colormap(ax12, turbo); caxis([0 1]); xticks([]); yticks([]); title(timePt); colorbar;

