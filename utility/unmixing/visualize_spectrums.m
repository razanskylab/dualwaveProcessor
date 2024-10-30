%% visualizes absorption spectrums of the following endogenous chromophores
% between 450nm and 1000nm:
% HbO2, HHb, lipid, melanin

% read HbO2 and HHb coefficients
mecHbFileId = fopen('C:\Code\dualwave-foam-dev\09_Utility\unmixing\absorptionSpectrum\hemoglobin.txt', 'r');
formatSpec = '%f';
mecHb = fscanf(mecHbFileId, formatSpec);
nRows = 3; % wavelength, HbO2, Hb
nCols = length(mecHb) / nRows;
mecHb = reshape(mecHb, [nRows, nCols]);
lambdas = mecHb(1,:); mec_HbO2 = mecHb(2,:); mec_HHb = mecHb(3,:);
% convert unit to [cm-1], see Prahl's page for details
mec_HbO2 = mec_HbO2 .* (150*2.303) ./ 64500;
mec_HHb = mec_HHb .* (150*2.303) ./ 64500;

% lipid
mecLipidFileId = fopen('C:\Code\dualwave-foam-dev\09_Utility\unmixing\absorptionSpectrum\lipid.txt', 'r');
formatSpec = '%f';
mecLipid = fscanf(mecLipidFileId, formatSpec);
nRows = 2; % wavelength, absorption coefficients
nCols = length(mecLipid) / nRows;
mecLipid = reshape(mecLipid, [nRows, nCols]);
mec_lipid = mecLipid(2,:);
mec_lipid = mec_lipid(1:2:end); % [m-1]
mec_lipid = mec_lipid ./ 1e2; % [cm-1]

% melanin
mec_melanin = 1.7e12 .* (lambdas.^(-3.48));

% plot the spectrums for all absorbers
figure; set(gcf, 'color', 'w');
plot(lambdas, log10(mec_HbO2), 'LineWidth', 2); hold on;
plot(lambdas, log10(mec_HHb), 'LineWidth', 2); hold on;
plot(lambdas, log10(mec_melanin), 'LineWidth', 2); hold on;
plot(lambdas, log10(mec_lipid), 'LineWidth', 2); hold off;
xlim([lambdas(1), lambdas(end)]); ylim([-3, 3]);
yticklabels({'0.001', '0.01', '0.1', '1', '10', '100', '1000'});
legend('oxygenated hemoglobin', 'deoxygenated hemoglobin', 'melanin', 'lipid');
grid on;
xlabel('wavelength [nm]');
ylabel('absorption coefficient [cm-1]');

% plot the spectrums for HbO2 and HHb
figure; set(gcf, 'color', 'w');
plot(lambdas, log10(mec_HbO2), 'LineWidth', 3, 'Color', [0.85, 0.325, 0.098]); hold on;
plot(lambdas, log10(mec_HHb), 'LineWidth', 3, 'Color', [0 0.447 0.741]); hold off;
xlim([500 590]); xticks([500 : 30 : 590]); 
ylim([1.8 2.5]); yticks([1.8 : 0.2 : 2.5]);
grid on;
legend('HbO2', 'HHb', 'FontSize', 16);
xlabel('wavelength [nm]', 'FontSize', 16);
ylabel('absorption coefficient [cm-1]', 'FontSize', 16);



