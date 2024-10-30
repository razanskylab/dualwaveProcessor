%% check the correspondence between power meter and photodiodes (pre-fiber and in-fiber)

data_base = 'C:\data_temp\2023_03_23';

currWavelength = 558;

if currWavelength == 532
    ppeList = [200, 500, 750, 900];
elseif currWavelength == 558
    ppeList = [200, 500, 700, 760];
else
    error('Wrong wavelength!');
end

prfList = [1, 2, 5, 7]; % [kHz]

ppeList = [500];
prfList = [5, 7];

nShotsAline = 50000; % in SNRScope, I acquired 50k shots at each PRF

currDir = fullfile(data_base, strcat('pdCalib', num2str(currWavelength)));

nPpe = length(ppeList);
nPrf = length(prfList);
% [x,y] coordinates for all shots per PRF, x is PD, y is PM 
shotPts_preFiber = zeros(nShotsAline*nPpe, 2, nPrf);
shotPts_inFiber = zeros(nShotsAline*nPpe, 2, nPrf);

for iPrf = 1 : nPrf
    currPrf = prfList(iPrf)
    
    for iPpe = 1 : nPpe
        currPpe = ppeList(iPpe);
        
        % read in corresponding power meter file
        currPMFileName = strcat('ppe', num2str(currPpe), '_prf', num2str(currPrf), 'k.csv');
        currPMFile = readmatrix(fullfile(currDir, currPMFileName));
        currPMMeas = currPMFile(:,2) .* 1e9; % [nJ], dim: [nShotsAline, 1]
        if length(currPMMeas) < nShotsAline
            error('Power meter did NOT get the right number of shots!');
        elseif length(currPMMeas) > nShotsAline
            error('Power meter got more shots than specified. How come?!');
        else
            disp('Power meter got just the right number of shots:)');
        end
        
        % read in corresponding photodiode file
        currPDFileName = strcat('ppe', num2str(currPpe), '_prf', num2str(currPrf), 'k.mat');
        currPDFile = load(fullfile(currDir, currPDFileName));
        preFiberMeas = squeeze(currPDFile.S.pdData);
        inFiberMeas = squeeze(currPDFile.S.usData);
        % correct laser jitter
        [preFiberMeas, ~, ~] = correct_laser_jitter(preFiberMeas, ones(size(preFiberMeas)));
        [inFiberMeas, ~, ~] = correct_laser_jitter(inFiberMeas, ones(size(inFiberMeas)));
        % convert counts to voltage
        preFiberMeas = counts_to_voltage(preFiberMeas, currPDFile.S.sensitivityPd, 16);
        inFiberMeas = counts_to_voltage(inFiberMeas, currPDFile.S.sensitivityUs, 16);

        if currPDFile.S.sensitivityPd ~= currPDFile.S.sensitivityUs
            error('Sensitivity difference detected!');
        end
        
        preFiberMeas = sum(preFiberMeas, 1);
        preFiberMeas = preFiberMeas(:); % dim: [nShotsAline, 1]
        inFiberMeas = sum(inFiberMeas, 1);
        inFiberMeas = inFiberMeas(:); % dim: [nShotsAline, 1]
        
        currIndices = (iPpe-1)*nShotsAline+1 : iPpe*nShotsAline;
        % assign to preFiber variable
        shotPts_preFiber(currIndices, 1, iPrf) = preFiberMeas;
        shotPts_preFiber(currIndices, 2, iPrf) = currPMMeas;
        % assign to inFiber variable
        shotPts_inFiber(currIndices, 1, iPrf) = inFiberMeas;
        shotPts_inFiber(currIndices, 2, iPrf) = currPMMeas;
        
    end
    
end

%% 

data_base = 'C:\data_temp\2023_05_02\pdCalib';

currWavelength = 558;
prfList = [2, 3, 4, 5, 6, 7]; % [kHz]
nPrf = length(prfList);
nShotsAline = 50000; % in SNRScope, I acquired 50k shots at each PRF
% [x,y] coordinates for all shots per PRF, x is PD, y is PM
shotPts_inFiber = zeros(nShotsAline, 2, nPrf);

for iPrf = 1 : nPrf
    currPrf = prfList(iPrf)
    
    if currWavelength == 558
        currPMFileName = strcat('dye', num2str(currWavelength), '_prf', num2str(currPrf), 'k.csv');
        currPDFileName = strcat('dye', num2str(currWavelength), '_prf', num2str(currPrf), 'k.mat');
    elseif currWavelength == 532
        currPMFileName = strcat('onda', num2str(currWavelength), '_prf', num2str(currPrf), 'k.csv');
        currPDFileName = strcat('onda', num2str(currWavelength), '_prf', num2str(currPrf), 'k.mat');
    else
        error('Invalid wavelength!');
    end
    
    % read in corresponding power meter file
    currPMFile = readmatrix(fullfile(data_base, currPMFileName));
    currPMMeas = currPMFile(:,2) .* 1e9; % [nJ], dim: [nShotsAline, 1]
    if length(currPMMeas) < nShotsAline
        error('Power meter did NOT get the right number of shots!');
    elseif length(currPMMeas) > nShotsAline
        error('Power meter got more shots than specified. How come?!');
    else
        disp('Power meter got just the right number of shots:)');
    end

    % read in corresponding photodiode file
    currPDFile = load(fullfile(data_base, currPDFileName));
    inFiberMeas = squeeze(currPDFile.S.pdData);
    % correct laser jitter
    [inFiberMeas, ~, ~] = correct_laser_jitter(inFiberMeas, ones(size(inFiberMeas)));
    % convert counts to voltage
    inFiberMeas = counts_to_voltage(inFiberMeas, currPDFile.S.sensitivityPd, 16);
    % estimate PPE in a.u.
    inFiberMeas = sum(inFiberMeas, 1);
    inFiberMeas = inFiberMeas(:); % dim: [nShotsAline, 1]
    % assign to shotPts variable
    shotPts_inFiber(:, 1, iPrf) = inFiberMeas;
    shotPts_inFiber(:, 2, iPrf) = currPMMeas;
    
end




%% generate scatter plots

% figure('Name', 'preFiber calibration plot');
% for iPrf = 1 : nPrf
%     scatter(squeeze(shotPts_preFiber(:,1,iPrf)), squeeze(shotPts_preFiber(:,2,iPrf)));
%     hold on;
% end
% 
% figure('Name', 'inFiber calibration plot');
% for iPrf = 1 : nPrf
%     scatter(squeeze(shotPts_inFiber(:,1,iPrf)), squeeze(shotPts_inFiber(:,2,iPrf)));
%     hold on;
% end

% % one plot for each PRF, color code photodiode
% for iPrf = 1 : nPrf
%     figure;
%     scatter(squeeze(shotPts_preFiber(:,1,iPrf)), squeeze(shotPts_preFiber(:,2,iPrf)), 'filled');
%     hold on;
%     scatter(squeeze(shotPts_inFiber(:,1,iPrf)), squeeze(shotPts_inFiber(:,2,iPrf)), 'filled');
%     hold off;
%     legend('preFiber PD', 'inFiber PD');
%     title(strcat('PRF: ', num2str(prfList(iPrf)), 'kHz'));
%     xlim([0 1200]); ylim([0 10e4]);
%     grid on;
% end


% tmpOrdered = tmp;
% tmpOrdered(2:2:end,:) = flip(tmpOrdered(2:2:end,:));
% tmpOrdered = tmpOrdered';
% tmpOrdered = tmpOrdered(:);