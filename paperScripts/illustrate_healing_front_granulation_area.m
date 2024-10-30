%% generate illustrative figure to show the area of:
% 1) wound healing front, i.e. band 1
% 2) granulation tissue area

clean_all;

data_base = 'D:\Roy\Experimental data\woundsO2Exp';
woundID = 'mouse1_r';
timePt = 'day4';

medfiltSize = 3;

pixelSize = 5e-6; % pixel size is 5 um
pixelArea = pixelSize * pixelSize;
woundDiam_day0 = 5e-3;
woundRadiIdx_day0 = round(woundDiam_day0/2/pixelSize);
nTheta = 1000;
deltaTheta = (2*pi) / nTheta;
theta = 0 : deltaTheta : 2*pi;

% Step 1: load dermis sO2 map as the base image
load(fullfile(data_base, timePt, 'processed\sO2', woundID));
nPixels = size(sO2_map_dermis, 1);

% Step 2: load generated band mask and granulation tissue mask
load(fullfile(data_base, timePt, 'processed\bandMasks', woundID));
load(fullfile(data_base, timePt, 'processed\granulationMask', woundID));

% Step 3: load generated wound boundary coordinates
currBoundary = ReadImageJROI(fullfile(data_base, timePt, 'processed\forMasking', strcat(woundID, '.roi')));
boundaryCoordinates = currBoundary.mnCoordinates;
geoCenter = round(mean(boundaryCoordinates, 1));
% calculate concentric boundaries by computing the vectors from each
% boundary coordinate and the geometric center, then push the boundary out
% in constant steps
diffVec = boundaryCoordinates - geoCenter;
diffVecLengths = vecnorm(diffVec,2,2);
diffVecDirections = diffVec ./ diffVecLengths;
bandWidth_mm = 0.5;
pixelSize_mm = 0.005;
bandWidth_pixels = round(bandWidth_mm/pixelSize_mm);
% band 1
bandIdx = 1;
band1Coordinates = round(geoCenter + diffVecDirections.*(diffVecLengths+bandWidth_pixels*bandIdx));
% take care of boundary condition at the image boundary
band1Coordinates(band1Coordinates<1) = 1;
band1Coordinates(band1Coordinates>nPixels) = nPixels;

xCoorCircle = cos(theta) .* woundRadiIdx_day0 + geoCenter(1);
yCoorCircle = sin(theta) .* woundRadiIdx_day0 + geoCenter(2);
% take care of boundary condition
xCoorCircle(xCoorCircle<1) = 1;
xCoorCircle(xCoorCircle>nPixels) = nPixels;
yCoorCircle(yCoorCircle<1) = 1;
yCoorCircle(yCoorCircle>nPixels) = nPixels;

figure; set(gcf, 'color', 'w');
imagesc(medfilt2(sO2_map_dermis'.*woundMask, [medfiltSize,medfiltSize]));
axis image; colormap turbo; caxis([0 1]);
xticks([]); yticks([]);
hold on;
scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 10, [1 1 1], 'filled'); 
hold on;
scatter(band1Coordinates(:,1), band1Coordinates(:,2), 10, [1 1 1], 'filled'); 
hold off;

figure; set(gcf, 'color', 'w');
imagesc(medfilt2(sO2_map_dermis'.*woundMask, [medfiltSize,medfiltSize]));
axis image; colormap turbo; caxis([0 1]);
xticks([]); yticks([]);
hold on;
scatter(boundaryCoordinates(:,1), boundaryCoordinates(:,2), 10, [1 1 1], 'filled'); 
hold on;
scatter(xCoorCircle, yCoorCircle, 10, [1,1,1], 'filled');
hold off;

% % write transparent band1 and granulation mask to disk
% alphachannel_band1 = ~band1Mask;
% imwrite(band1Mask, 'C:\Users\johan\Desktop\band1MaskTransparent_day10.png', 'Alpha', double(alphachannel_band1));
% 
% alphachannel_granu = ~granulationMask;
% imwrite(granulationMask, 'C:\Users\johan\Desktop\granuMaskTransparent_day10.png', 'Alpha', double(alphachannel_granu));




