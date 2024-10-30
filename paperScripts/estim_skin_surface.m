%% Script to estimate skin surface after wrap signal cropping
clean_all;

data_base = 'D:\Roy\Experimental data\woundsO2Exp\day-3';
data_name = 'mouse1_l_10um_us.mat';

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

%% Load pre-computed wrap cropping depths and apply onto US volume
mouseInfo = data_name(1:8);
load(fullfile(data_base, 'wrapCropDepths', mouseInfo));

usVol = zeros(size(dwp.usVol_pe));
[~, nAlines, nBscan] = size(usVol);
bscanInterval_mm = 0.25;
stepSize_mm = dwp.scanSett.dr(2);
bscanInterval_idx = bscanInterval_mm / stepSize_mm;
zVec_mm = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
zVec_mm = zVec_mm(dwp.zCropIdx_us(1):dwp.zCropIdx_us(2));
xVec_mm = (0:nAlines-1) .* stepSize_mm;
for iBscan = 1 : nBscan
    for iAline = 1 : nAlines
        
        currStartDepth_mm = wrapCropDepths_mm(iAline, iBscan);
        [~, currStatrtDepth_idx] = min(abs(zVec_mm - currStartDepth_mm));
        usVol(currStatrtDepth_idx:end,iAline,iBscan) = dwp.usVol_pe(currStatrtDepth_idx:end,iAline,iBscan);
        
    end
end
% visualize B-scans of wrap-cropped US volume
fBscan = figure('Name', 'B-scan', 'Color', 'w');
fBscan.Units = 'normalized';
fBscan.OuterPosition = [0 0 1 1];
for iBscan = 1 : bscanInterval_idx : nBscan
    currBscan = squeeze(usVol(:,:,iBscan));
    imagesc(xVec_mm, zVec_mm, currBscan);
    colormap redblue; colorbar; axis tight equal;
    caxis([-max(abs(currBscan(:))), max(abs(currBscan(:)))]);
    caxis([-100 100]);
    xlabel('[mm]'); ylabel('[mm]');
    title(sprintf('%i-th B-scan', iBscan));
    waitforbuttonpress;
end

%% estimate skin surface based on wrap-cropped US volume
noiseAmplitude = 4.25;
smoothRange_mm = 1;
medfiltWindowSize = smoothRange_mm / stepSize_mm;

noiseAmplitudeAdjusted = noiseAmplitude - 0;
noiseAmplitudeAdjustAlineIdx = 401;
noiseAmplitudeAdjustBscanIdx = 401;

skinStartIndices_bscan = zeros(1, nAlines);
skinStartDepths_mm_bscan = zeros(1, nAlines);
% allocate memory for crop indices of the entire FOV
skinCropDepths_mm = zeros(nAlines, nBscan);

for iBscan = 1 : nBscan
    currBscan = usVol(:, :, iBscan);
    currBscan = medfilt2(currBscan, [3,3]);
    
    for iAline = 1 : nAlines
        currAline = currBscan(:, iAline);
        
        if (iAline >= noiseAmplitudeAdjustAlineIdx) & (iBscan>=noiseAmplitudeAdjustBscanIdx)
            sigIndices = find(currAline < -noiseAmplitudeAdjusted);
        else
            sigIndices = find(currAline < -noiseAmplitude);
        end
        
        if isempty(sigIndices)
            skinStartIndices_bscan(iAline) = nan;
            skinStartDepths_mm_bscan(iAline) = nan;
        else
            skinStartIdx = sigIndices(1);
            skinStartDepth_mm = zVec_mm(skinStartIdx);
            skinStartIndices_bscan(iAline) = skinStartIdx;
            skinStartDepths_mm_bscan(iAline) = skinStartDepth_mm;
        end
        
    end
    
    smoothedSkinProfile = movmedian(skinStartDepths_mm_bscan, medfiltWindowSize, 'omitnan');
    skinCropDepths_mm(:,iBscan) = smoothedSkinProfile;
    
    % show some B-scans together with the cropping profile
    if ~mod(iBscan,bscanInterval_idx)
        fBscan = figure('Name', 'B-scan', 'Color', 'w');
        fBscan.Units = 'normalized';
        fBscan.OuterPosition = [0 0 1 1];
        imagesc(xVec_mm, zVec_mm, currBscan);
        colormap redblue; colorbar; axis tight equal;
        caxis([-max(abs(currBscan(:))), max(abs(currBscan(:)))]);
        xticks([]); yticks([]);
        xlabel('[mm]', 'FontSize', 14); ylabel('[mm]', 'FontSize', 14);
        title(sprintf('%i-th B-scan', iBscan));
        hold on;
        scatter(xVec_mm, smoothedSkinProfile, 'g', 'filled');
        waitforbuttonpress;
        close;
    end
    
end

% save skin surface depths
save_name = data_name(1:8); % selects mouse number and left/right
save_dir = fullfile(data_base, 'skinCropDepths');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
% save(fullfile(save_dir, save_name), 'skinCropDepths_mm');

%% visualize estimated skin surface
close all;

figure; imagesc(skinCropDepths_mm); axis image; colorbar;
title('Check skin surface smoothness');

fSurf = figure('Name', 'Surface plot', 'Color', 'w');
fSurf.Units = 'normalized';
fSurf.OuterPosition = [0 0 1 1];
yVec_mm = (0:nBscan-1) .* stepSize_mm;
surf(yVec_mm, xVec_mm, skinCropDepths_mm);
set(gca, 'Zdir', 'reverse');
xlabel('y [mm]', 'FontSize',16); ylabel('x [mm]', 'FontSize',16); zlabel('z [mm]', 'FontSize',16);
title('Estimated skin surface', 'FontSize',16);

% surface mesh plot
[X,Y] = meshgrid(0:700);
skinSurfaceMesh = skinCropDepths_mm;
figure; mesh(X,Y,skinSurfaceMesh); 
set(gca, 'Zdir', 'reverse');
caxis([6.8 7.8]); colorbar;
xticks([]); yticks([]); set(gca, 'Visible', 'Off');
colormap bone;

% surface mesh plot - trying downsampling
skinSurfaceMesh = skinCropDepths_mm(1:7:701, 1:7:701);
[X,Y] = meshgrid(0:100);
figure; mesh(X,Y,skinSurfaceMesh); 
set(gca, 'Zdir', 'reverse');
caxis([6.8 7.8]); colorbar;
xticks([]); yticks([]); set(gca, 'Visible', 'Off');
colormap bone;

