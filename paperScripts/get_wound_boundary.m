%% script to delineate wound boundary and generate concentric boundaries around it
clean_all;

data_base = 'D:\Roy\Experimental data\woundsO2Exp';
woundID = 'mouse1_r';
timePt = 'day4';
load(fullfile(data_base, timePt, 'processed/structure', woundID));
nPixels = size(structure_map_allDepths,1);

% export native picture to paint boundary in Fiji
figure; imagesc(structure_map_allDepths'); axis image;
colormap bone; caxis([0 80]); xticks([]); yticks([]); hold on;
export_fig(fullfile(data_base, timePt, 'processed\forMasking', woundID), '-native', '-tiff');

% check if the work in Fiji is done
prompt = "Have you saved the hand-drawn ROI from Fiji?";
reply = input(prompt, "s");
if strcmp(reply, "yes")
    disp('Yay! We can proceed.');
else
    error('Please draw wound boundary in Fiji and save the ROI file!');
end

% verify if the painted boundary is good enough (covers wound area and convex)
currBoundary = ReadImageJROI(fullfile(data_base, timePt, 'processed\forMasking', strcat(woundID, '.roi')));
boundaryCoordinates = currBoundary.mnCoordinates;
scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 5, 'g', 'filled'); hold on;
prompt = "Are you happy with the hand-drawn ROI from Fiji?";
reply = input(prompt, "s");
if strcmp(reply, "yes")
    disp('Yay! We can proceed.');
else
    error('Please re-draw wound boundary in Fiji and save the ROI file!');
end

% once we are happy with the painted wound boundary, proceed with
% generating concentric bands around the boundary
geoCenter = round(mean(boundaryCoordinates, 1));
scatter(geoCenter(1), geoCenter(2), 'r', 'filled'); hold on;

% calculate concentric boundaries by computing the vectors from each
% boundary coordinate and the geometric center, then push the boundary out
% in constant steps
diffVec = boundaryCoordinates - geoCenter;
diffVecLengths = vecnorm(diffVec,2,2);
diffVecDirections = diffVec ./ diffVecLengths;
bandWidth_mm = 1;
pixelSize_mm = 0.005;
bandWidth_pixels = round(bandWidth_mm/pixelSize_mm);

% band 1
bandIdx = 1;
band1Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
% take care of boundary condition at the image boundary
band1Coordinates(band1Coordinates<1) = 1;
band1Coordinates(band1Coordinates>nPixels) = nPixels;
scatter(band1Coordinates(:,1), band1Coordinates(:,2), 5, [0.9290 0.6940 0.1250], 'filled'); hold on;

% band 2
bandIdx = 2;
band2Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
% take care of boundary condition at the image boundary
band2Coordinates(band2Coordinates<1) = 1;
band2Coordinates(band2Coordinates>nPixels) = nPixels;
scatter(band2Coordinates(:,1), band2Coordinates(:,2), 5, [0.8500 0.3250 0.0980], 'filled'); hold on;

%% get boundary and band masks
% once the boundary is verified in the overlay plot, proceed with
% generating binary masks
[xq, yq] = meshgrid(1:nPixels, 1:nPixels);
xq = xq(:); yq = yq(:);
woundMask = inpolygon(xq, yq, boundaryCoordinates(:,1), boundaryCoordinates(:,2));
woundMask = reshape(~woundMask,[nPixels,nPixels]);
figure; imagesc(woundMask); axis image; h = colorbar; h.Ticks = [0, 1];
xticks([]); yticks([]);

band1Mask = inpolygon(xq, yq, band1Coordinates(:,1), band1Coordinates(:,2));
band1Mask = reshape(band1Mask,[nPixels,nPixels]);
band1Mask = band1Mask .* woundMask;
figure; imagesc(band1Mask); axis image; h = colorbar; h.Ticks = [0, 1];
xticks([]); yticks([]);

band2Mask = inpolygon(xq, yq, band2Coordinates(:,1), band2Coordinates(:,2));
band2Mask = reshape(band2Mask,[nPixels,nPixels]);
band2Mask = band2Mask .* woundMask .* (~band1Mask);
figure; imagesc(band2Mask); axis image; h = colorbar; h.Ticks = [0, 1];
xticks([]); yticks([]);

remainingMask = woundMask .* (~band1Mask) .* (~band2Mask);
figure; imagesc(remainingMask); axis image; h = colorbar; h.Ticks = [0, 1];
xticks([]); yticks([]);

%% make composite structural and sO2 image

load(fullfile(data_base, timePt, 'processed/sO2', woundID));

medfiltSize = 3;

% structural composite
structure_map_composite = (medfilt2(structure_map_dermis',[medfiltSize,medfiltSize]).*band1Mask) + ...
                          (medfilt2(structure_map_hypodermis',[medfiltSize,medfiltSize]).*band2Mask) + ...
                          (medfilt2(structure_map_allDepths',[medfiltSize,medfiltSize]).*remainingMask);
figure; imagesc(structure_map_composite); axis image;
colormap bone; caxis([0 50]); xticks([]); yticks([]); hold on;
scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 5, 'g', 'filled'); hold on;
scatter(geoCenter(1), geoCenter(2), 'r', 'filled'); hold on;
scatter(band1Coordinates(:,1), band1Coordinates(:,2), 5, [0.9290 0.6940 0.1250], 'filled'); hold on;
scatter(band2Coordinates(:,1), band2Coordinates(:,2), 5, [0.8500 0.3250 0.0980], 'filled'); hold off;
                    
% sO2 composite
sO2_map_composite = (medfilt2(sO2_map_dermis',[medfiltSize,medfiltSize]).*band1Mask) + ...
                    (medfilt2(sO2_map_hypodermis',[medfiltSize,medfiltSize]).*band2Mask) + ...
                    (medfilt2(sO2_map_allDepths',[medfiltSize,medfiltSize]).*remainingMask);
figure; imagesc(sO2_map_composite); axis image;
colormap turbo; caxis([0 1]); xticks([]); yticks([]);

%% Visualize long-term sO2 dynamics with wound masked out

data_base = "C:\data_temp\woundExp";
woundID = "mouse3_l";
timePt_list = ["day-3", "day4", "day7", "day10"];

nPixels = 1401;
[xq, yq] = meshgrid(1:nPixels, 1:nPixels);
xq = xq(:); yq = yq(:);

figure('Name', woundID); set(gcf, 'color', 'w');
medfiltSize = 3;

nTimePt = length(timePt_list);
for iTimePt = 1 : nTimePt
    currTimePt = timePt_list(iTimePt)
    
    currSO2 = load(fullfile(data_base, currTimePt, "processed\sO2", woundID));
    
    if iTimePt > 1
        % just for mouse3_l cause US dataset is missing at day4
        if iTimePt == 2
            continue;
        end
        currBoundary = ReadImageJROI(fullfile(data_base, currTimePt, "processed\forMasking", strcat(woundID, '.roi')));
        boundaryCoordinates = currBoundary.mnCoordinates;
        woundMask = inpolygon(xq, yq, boundaryCoordinates(:,1), boundaryCoordinates(:,2));
        woundMask = reshape(~woundMask,[nPixels,nPixels]);
        sO2_map_dermis = currSO2.sO2_map_dermis' .* woundMask;
        sO2_map_dermis = medfilt2(sO2_map_dermis, [medfiltSize,medfiltSize]);
        sO2_map_hypodermis = currSO2.sO2_map_hypodermis' .* woundMask;
        sO2_map_hypodermis = medfilt2(sO2_map_hypodermis, [medfiltSize,medfiltSize]);
        subplot(2,nTimePt,iTimePt); imagesc(sO2_map_dermis);
        axis image; xticks([]); yticks([]); colormap turbo; caxis([0 1]); title(currTimePt);
        hold on; scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 3, 'w', 'filled'); hold off;
        subplot(2,nTimePt,nTimePt+iTimePt); imagesc(sO2_map_hypodermis);
        axis image; xticks([]); yticks([]); colormap turbo; caxis([0 1]); title(currTimePt);
        hold on; scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 3, 'w', 'filled'); hold off;
    else
        subplot(2,nTimePt,iTimePt); imagesc(medfilt2(currSO2.sO2_map_dermis',[medfiltSize,medfiltSize]));
        axis image; xticks([]); yticks([]); colormap turbo; caxis([0 1]); title(currTimePt);
        subplot(2,nTimePt,nTimePt+iTimePt); imagesc(medfilt2(currSO2.sO2_map_hypodermis',[medfiltSize,medfiltSize]));
        axis image; xticks([]); yticks([]); colormap turbo; caxis([0 1]); title(currTimePt);
    end
    
end

%% generate masks for all mice and time points
clean_all;

bandWidth_mm = 0.5;
pixelSize_mm = 0.005;
nBands = 6; % 0, 0.5, 1, 1.5, 2mm and remaining pixels

% read the dermis structural and sO2 image
data_base = 'C:\data_temp\woundExp';
woundID_list = ["mouse1_l", "mouse1_r", "mouse2_r", "mouse3_r", "mouse4_l", "mouse4_r"];
timePt_list = ["day4", "day7", "day10"];

woundID_list = ["mouse1_r"];
timePt_list = ["day4"];

for iWound = 1 : length(woundID_list)
    woundID = woundID_list(iWound)
    
    for iTimePt = 1 : length(timePt_list)
        timePt = timePt_list(iTimePt)

        structure = load(fullfile(data_base, timePt, 'processed/structure', woundID));
        sO2 = load(fullfile(data_base, timePt, 'processed/sO2', woundID));

        figure; set(gcf, 'color', 'w');
        ax1 = subplot(121); imagesc(sO2.sO2_map_dermis'); axis image;
        colormap(ax1, turbo); caxis([0 1]);
        ax2 = subplot(122); imagesc(structure.structure_map_dermis'); axis image;
        colormap(ax2, bone); caxis([0 50]); hold on;

        nPixels = size(structure.structure_map_allDepths,1);
        [xq, yq] = meshgrid(1:nPixels, 1:nPixels);
        xq = xq(:); yq = yq(:);

        % load the saved ROI file
        currBoundary = ReadImageJROI(fullfile(data_base, timePt, 'processed\forMasking', strcat(woundID, '.roi')));
        boundaryCoordinates = currBoundary.mnCoordinates;

        % check the wound boundary by overlaying it on the structure and sO2 image
        scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 3, 'g', 'filled'); hold off;

        % generate band masks

        % band 0: wound mask
        bandIdx = 0;
        woundMask = inpolygon(xq, yq, boundaryCoordinates(:,1), boundaryCoordinates(:,2));
        woundMask = reshape(~woundMask,[nPixels,nPixels]);

        figure; set(gcf, 'color', 'w');
        subplot(1, nBands, bandIdx+1); imagesc(woundMask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('wound mask');

        geoCenter = round(mean(boundaryCoordinates, 1));
        diffVec = boundaryCoordinates - geoCenter;
        diffVecLengths = vecnorm(diffVec,2,2);
        diffVecDirections = diffVec ./ diffVecLengths;
        bandWidth_pixels = round(bandWidth_mm/pixelSize_mm);

        % band 1
        bandIdx = 1;
        band1Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
        % take care of boundary condition at the image boundary
        band1Coordinates(band1Coordinates<1) = 1;
        band1Coordinates(band1Coordinates>nPixels) = nPixels;
        band1Mask = inpolygon(xq, yq, band1Coordinates(:,1), band1Coordinates(:,2));
        band1Mask = reshape(band1Mask,[nPixels,nPixels]);
        band1Mask = band1Mask .* woundMask;
        subplot(1, nBands, bandIdx+1); imagesc(band1Mask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('band 1 mask'); % 0.5mm out of wound boundary

        % band 2
        bandIdx = 2;
        band2Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
        % take care of boundary condition at the image boundary
        band2Coordinates(band2Coordinates<1) = 1;
        band2Coordinates(band2Coordinates>nPixels) = nPixels;
        band2Mask = inpolygon(xq, yq, band2Coordinates(:,1), band2Coordinates(:,2));
        band2Mask = reshape(band2Mask,[nPixels,nPixels]);
        band2Mask = band2Mask .* woundMask .* (~band1Mask);
        subplot(1, nBands, bandIdx+1); imagesc(band2Mask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('band 2 mask'); % 1mm out of wound boundary

        % band 3
        bandIdx = 3;
        band3Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
        % take care of boundary condition at the image boundary
        band3Coordinates(band3Coordinates<1) = 1;
        band3Coordinates(band3Coordinates>nPixels) = nPixels;
        band3Mask = inpolygon(xq, yq, band3Coordinates(:,1), band3Coordinates(:,2));
        band3Mask = reshape(band3Mask,[nPixels,nPixels]);
        band3Mask = band3Mask .* woundMask .* (~band1Mask) .* (~band2Mask);
        subplot(1, nBands, bandIdx+1); imagesc(band3Mask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('band 3 mask'); % 1.5mm out of wound boundary

        % band 4
        bandIdx = 4;
        band4Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
        % take care of boundary condition at the image boundary
        band4Coordinates(band4Coordinates<1) = 1;
        band4Coordinates(band4Coordinates>nPixels) = nPixels;
        band4Mask = inpolygon(xq, yq, band4Coordinates(:,1), band4Coordinates(:,2));
        band4Mask = reshape(band4Mask,[nPixels,nPixels]);
        band4Mask = band4Mask .* woundMask .* (~band1Mask) .* (~band2Mask) .* (~band3Mask);
        subplot(1, nBands, bandIdx+1); imagesc(band4Mask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('band 4 mask'); % 2mm out of wound boundary

        % band 5 - remaining pixels
        bandIdx = 5;
        band5Mask = woundMask .* (~band1Mask) .* (~band2Mask) .* (~band3Mask) .* (~band4Mask);
        subplot(1, nBands, bandIdx+1); imagesc(band5Mask); axis image; h = colorbar; h.Ticks = [0, 1];
        xticks([]); yticks([]); title('band 5 mask'); % remaining pixels

        % save the generated masks
        save(fullfile(data_base, timePt, 'processed\bandMasks', woundID), ...
                'woundMask', 'band1Mask', 'band2Mask', 'band3Mask', 'band4Mask', 'band5Mask');

    end
end

