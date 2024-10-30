%% Script to screen depth of wrap signals and crop them away based on US dataset
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

%% Check a few B-scans to have an idea on the depth of wrap in the entire FOV
nBscan = size(dwp.usVol_pe,3);
bscanInterval_mm = 0.25;
stepSize_mm = dwp.scanSett.dr(2);
bscanInterval_idx = bscanInterval_mm / stepSize_mm;

zVec_mm = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
zVec_mm = zVec_mm(dwp.zCropIdx_us(1):dwp.zCropIdx_us(2));
nAlines = size(dwp.usVol_pe,2);
xVec_mm = (0:nAlines-1) .* stepSize_mm;

fBscan = figure('Name', 'B-scan', 'Color', 'w');
fBscan.Units = 'normalized';
fBscan.OuterPosition = [0 0 1 1];
for iBscan = 1 : bscanInterval_idx : nBscan
    currBscan = squeeze(dwp.usVol_pe(:,:,iBscan));
    imagesc(xVec_mm, zVec_mm, currBscan);
    colormap redblue; colorbar; axis tight equal;
%     caxis([-max(abs(currBscan(:))), max(abs(currBscan(:)))]);
    caxis([-100 100]);
    xlabel('[mm]'); ylabel('[mm]');
    title(sprintf('%i-th B-scan', iBscan));
    waitforbuttonpress;
end

%% If there's an easy linear depth cropping, we are lucky
wrapCropDepths_mm = zeros(nAlines, nBscan) + 6.6;
% visualize linear cut
for iBscan = 1 : nBscan
    currBscan = dwp.usVol_pe(:, :, iBscan);
    currBscan = medfilt2(currBscan, [3,3]);

    % show some B-scans together with the cropping profile
    if (~mod(iBscan,bscanInterval_idx) || (iBscan==1))
        fBscan = figure('Name', 'B-scan', 'Color', 'w');
        fBscan.Units = 'normalized';
        fBscan.OuterPosition = [0 0 1 1];
        imagesc(xVec_mm, zVec_mm, currBscan);
        colormap redblue; colorbar; axis tight equal;
%         caxis([-max(abs(currBscan(:))), max(abs(currBscan(:)))]);
        caxis([-100 100]);
        xlabel('[mm]'); ylabel('[mm]');
        title(sprintf('%i-th B-scan', iBscan));
        hold on;
        linearCrop_mm = wrapCropDepths_mm(:,iBscan);
        scatter(xVec_mm, linearCrop_mm, 'g', 'filled');
        waitforbuttonpress;
        close;
    end
    
end

% save wrap cropping surface depths
save_name = data_name(1:8); % selects mouse number and left/right
save_dir = fullfile(data_base, 'wrapCropDepths');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save(fullfile(save_dir, save_name), 'wrapCropDepths_mm');

%% If there's NO easy linear cropping, estimate the wrap profile and crop based on that
%  detect wrap signal depth per B-scan
skinHighestPt_mm = 6.75; % means any signal above this depth is wrap
noiseAmplitude = 5.5; % above this number is considered legit US signal
smoothRange_mm = 0.25; % assumes the wrap surface is flat within this range
movDownSize_mm = 0.042; % move down this much to crop away wrap signal
bscanSeparationIdx = 1000; % may need to process chuncks of B-scans separately, depending on each dataset
goodxMarginEnd_mm = 3; % along one b-scan the separation between wrap and skin can also be very different
goodxMarginEnd_idx = goodxMarginEnd_mm / stepSize_mm;
medfiltWindowSize = smoothRange_mm / stepSize_mm;
% fillmissingWindow_mm = 1;
% fillMissingWindowSize = fillmissingWindow_mm / stepSize_mm;

wrapStartIndices_bscan = zeros(1, nAlines);
wrapStartDepths_mm_bscan = zeros(1, nAlines);

% allocate memory for crop indices of the entire FOV
wrapCropDepths_mm = zeros(nAlines, nBscan);

for iBscan = 1 : nBscan
    currBscan = dwp.usVol_pe(:, :, iBscan);
    currBscan = medfilt2(currBscan, [3,3]);
    
    for iAline = 1 : nAlines
        
        currAline = currBscan(:, iAline);
%         sigIndices = find(currAline < -noiseAmplitude);
        sigIndices = find(abs(currAline) > noiseAmplitude);
        if isempty(sigIndices)
            wrapStartIndices_bscan(iAline) = nan;
            wrapStartDepths_mm_bscan(iAline) = nan;
        else
            wrapStartIdx = sigIndices(1);
            wrapStartDepth_mm = zVec_mm(wrapStartIdx);
            wrapStartIndices_bscan(iAline) = wrapStartIdx;
            wrapStartDepths_mm_bscan(iAline) = wrapStartDepth_mm;
        end
        
    end
    
    smoothedWrapProfile = movmedian(wrapStartDepths_mm_bscan, medfiltWindowSize, 'omitnan');
    smoothedWrapProfile = smoothedWrapProfile + movDownSize_mm;
    % get rid of falsed detection of actual skin surface due to low
    % wrap signal for B-scans defined by bscanSeparationIdx (after this index)
    if iBscan <= bscanSeparationIdx
        smoothedWrapProfileCorrected = smoothedWrapProfile;
    else
        smoothedWrapProfileCorrected = smoothedWrapProfile;
%         smoothedWrapProfileCorrected(smoothedWrapProfileCorrected>skinHighestPt_mm) = nan;
%         smoothedWrapProfileCorrected(goodxMarginEnd_idx+1:end) = smoothedWrapProfile(goodxMarginEnd_idx+1:end);
%         fillConst = mean(smoothedWrapProfileCorrected, 'omitnan');
%         smoothedWrapProfileCorrected = fillmissing(smoothedWrapProfileCorrected, 'constant', fillConst);
        smoothedWrapProfileCorrected(goodxMarginEnd_idx+1:end) = nan;
%         fillConst = mean(smoothedWrapProfileCorrected, 'omitnan');
        fillConst = smoothedWrapProfileCorrected(goodxMarginEnd_idx);
        smoothedWrapProfileCorrected = fillmissing(smoothedWrapProfileCorrected, 'constant', fillConst);
    end
    
    wrapCropDepths_mm(:,iBscan) = smoothedWrapProfileCorrected;
    
    % show some B-scans together with the cropping profile
    if ~mod(iBscan,bscanInterval_idx)
        fBscan = figure('Name', 'B-scan', 'Color', 'w');
        fBscan.Units = 'normalized';
        fBscan.OuterPosition = [0 0 1 1];
        imagesc(xVec_mm, zVec_mm, currBscan);
        colormap redblue; colorbar; axis tight equal;
        caxis([-max(abs(currBscan(:))), max(abs(currBscan(:)))]);
        xlabel('[mm]'); ylabel('[mm]');
        title(sprintf('%i-th B-scan', iBscan));
        hold on;
%         scatter(xVec_mm, smoothedWrapProfileCorrected, 'g', 'filled');
        plot(xVec_mm, smoothedWrapProfileCorrected, 'g', 'LineWidth', 2);
        waitforbuttonpress;
        close;
    end
    
end

% save wrap cropping surface depths
save_name = data_name(1:8); % selects mouse number and left/right
save_dir = fullfile(data_base, 'wrapCropDepths');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save(fullfile(save_dir, save_name), 'wrapCropDepths_mm');

%% visualize wrap crop surface
close all;

figure; imagesc(wrapCropDepths_mm); axis image; colorbar;
title('Check wrap crop surface smoothness');

fSurf = figure('Name', 'Surface plot', 'Color', 'w');
fSurf.Units = 'normalized';
fSurf.OuterPosition = [0 0 1 1];
yVec_mm = (0:nBscan-1) .* stepSize_mm;
surf(yVec_mm, xVec_mm, wrapCropDepths_mm);
set(gca, 'Zdir', 'reverse');
xlabel('y [mm]', 'FontSize',16); ylabel('x [mm]', 'FontSize',16); zlabel('z [mm]', 'FontSize',16);
title('Surface to crop out wrap signal', 'FontSize',16);



