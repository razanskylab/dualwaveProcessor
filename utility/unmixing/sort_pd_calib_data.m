% Sort out Photodiode calibration data
% Power meter data is in each .csv file
% Photodiode readings are stored in the "SNRScope" object

% For each wavelength:
% put PD readings at all energy levels in matrix X, each column corresponds to each energy level
% put power meter readings at all energy levels in matrix Y, each column corresponds to each energy level

% clear, clc; 
% close all;

data_base = 'C:\data_temp\2022_08_04';

dye_pm_fileInfo = dir(fullfile(data_base, 'dye_*.csv'));
onda532_pm_fileInfo = dir(fullfile(data_base, 'onda532_*.csv'));

nEnergyLevels = size(dye_pm_fileInfo, 1);
assert(size(onda532_pm_fileInfo, 1) == nEnergyLevels);

nShots = 2000;

pmMat_dye = zeros(nShots, nEnergyLevels); % [nJ]
pmMat_onda532 = zeros(nShots, nEnergyLevels); % [nJ]

for i = 1 : nEnergyLevels

  currName_dye_pm = dye_pm_fileInfo(i).name;
  currMat_dye_pm = readmatrix(fullfile(data_base, currName_dye_pm));
  pmMat_dye(:,i) = currMat_dye_pm(:,2) .* 1e9;

  currName_onda532_pm = onda532_pm_fileInfo(i).name;
  currMat_onda532_pm = readmatrix(fullfile(data_base, currName_onda532_pm));
  pmMat_onda532(:,i) = currMat_onda532_pm(:,2) .* 1e9;

end

save(fullfile(data_base, 'pmData'), 'pmMat_dye', 'pmMat_onda532');


%% Process and save PD data

dye_pd_fileInfo = dir(fullfile(data_base, 'dye_*.mat'));
onda532_pd_fileInfo = dir(fullfile(data_base, 'onda532_*.mat'));

pdMat_dye = zeros(nShots, nEnergyLevels); % [a.u.]
pdMat_onda532 = zeros(nShots, nEnergyLevels); % [a.u.]

for i = 1 : nEnergyLevels

  currName_dye_pd = dye_pd_fileInfo(i).name;
  currS_dye_pd = load(fullfile(data_base, currName_dye_pd));
  pdShots_dye = squeeze(currS_dye_pd.S.pdData);
  usShots_dye = squeeze(currS_dye_pd.S.usData);
  [pdShots_dye, usShots_dye, ~] = correct_laser_jitter(pdShots_dye, usShots_dye);
  pdShots_dye = counts_to_voltage(pdShots_dye, currS_dye_pd.S.sensitivityPd, 16);
  pdShots_dye = pd_shot_to_value(pdShots_dye, [1,100]);
  pdMat_dye(:,i) = pdShots_dye(:);

  currName_onda532_pd = onda532_pd_fileInfo(i).name;
  currS_onda532_pd = load(fullfile(data_base, currName_onda532_pd));
  pdShots_onda532 = squeeze(currS_onda532_pd.S.pdData);
  usShots_onda532 = squeeze(currS_onda532_pd.S.usData);
  [pdShots_onda532, usShots_onda532, ~] = correct_laser_jitter(pdShots_onda532, usShots_onda532);
  pdShots_onda532 = counts_to_voltage(pdShots_onda532, currS_onda532_pd.S.sensitivityPd, 16);
  pdShots_onda532 = pd_shot_to_value(pdShots_onda532, [1,100]);
  pdMat_onda532(:,i) = pdShots_onda532(:);

end

save(fullfile(data_base, 'pdData'), 'pdMat_dye', 'pdMat_onda532');


%% compute correlation coefficient between PD and PM data
[rho_dye, ~] = corr(pdMat_dye, pmMat_dye);
[rho_onda532, ~] = corr(pdMat_onda532, pmMat_onda532);

% diag(rho_dye)
% diag(rho_onda532)

%% generate linear fit plot (PD as x-axis, PM as y-axis)

colorCode_dye = wavelength2color(558); % good Urs legacy:)
colorCode_onda532 = wavelength2color(532);

fFit = figure('Name', 'PD against PM fit', 'Color', 'w');
fFit.Units = 'normalized';
fFit.OuterPosition = [0 0 1 1];

pd_dye_median = median(pdMat_dye, 1);
pd_dye_LB = prctile(pdMat_dye, 25, 1);
pd_dye_UB = prctile(pdMat_dye, 75, 1);
pd_dye_neg = pd_dye_median - pd_dye_LB;
pd_dye_pos = pd_dye_UB - pd_dye_median;

pm_dye_median = median(pmMat_dye, 1);
pm_dye_LB = prctile(pmMat_dye, 25, 1);
pm_dye_UB = prctile(pmMat_dye, 75, 1);
pm_dye_neg = pm_dye_median - pm_dye_LB;
pm_dye_pos = pm_dye_UB - pm_dye_median;

errorbar(pd_dye_median, pm_dye_median, ...
         pm_dye_neg, pm_dye_pos, pd_dye_neg, pd_dye_pos, ...
         'LineWidth', 3, 'Color', colorCode_dye);

hold on;

pd_onda532_median = median(pdMat_onda532, 1);
pd_onda532_LB = prctile(pdMat_onda532, 25, 1);
pd_onda532_UB = prctile(pdMat_onda532, 75, 1);
pd_onda532_neg = pd_onda532_median - pd_onda532_LB;
pd_onda532_pos = pd_onda532_UB - pd_onda532_median;

pm_onda532_median = median(pmMat_onda532, 1);
pm_onda532_LB = prctile(pmMat_onda532, 25, 1);
pm_onda532_UB = prctile(pmMat_onda532, 75, 1);
pm_onda532_neg = pm_onda532_median - pm_onda532_LB;
pm_onda532_pos = pm_onda532_UB - pm_onda532_median;

errorbar(pd_onda532_median, pm_onda532_median, ...
         pm_onda532_neg, pm_onda532_pos, pd_onda532_neg, pd_onda532_pos, ...
         'LineWidth', 3, 'Color', colorCode_onda532);

hold on;

%% fit linear curve

fitOrder = 1;
polyFit_dye = polyfit(pd_dye_median, pm_dye_median, fitOrder);
polyFit_onda532 = polyfit(pd_onda532_median, pm_onda532_median, fitOrder);

% [polyFit_dye, errStruct_dye, ~] = polyfit(pd_dye_median, pm_dye_median, fitOrder);
% [polyFit_onda532, errStruct_onda532, ~] = polyfit(pd_onda532_median, pm_onda532_median, fitOrder);
% errStruct_dye.normr
% errStruct_onda532.normr

fitCurve_dye = polyval(polyFit_dye, pd_dye_median);
fitCurve_onda532 = polyval(polyFit_onda532, pd_onda532_median);

plot(pd_dye_median, fitCurve_dye, 'm--', 'LineWidth', 3);

hold on;

plot(pd_onda532_median, fitCurve_onda532, 'm--', 'LineWidth', 3);

hold off;

grid on;
xlabel('Photodiode readings [a.u.]', 'FontSize', 16);
ylabel('Power meter readings [nJ]', 'FontSize', 16);
xlim([4000 12000]); ylim([450 950]);
title('PD calibration plot', 'FontSize', 16);
legend('Dye - measured', 'Onda532 - measured', ...
       'Dye - fitted', 'Onda532 - fitted', 'FontSize', 16);

