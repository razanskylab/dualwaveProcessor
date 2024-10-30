%% compute the determinant of molar extinction coefficient matrices
%  one fixed wavelength: 532 nm
%  the other dye wavelength range: 534 - 610 nm

proj_base = 'C:\Code\dualwave-processor';
extinction_file = 'Molar_Extinction_Coefficients_Hemoglobin_Water_532_to_610.txt';
mecFileId = fopen(fullfile(proj_base, extinction_file), 'r');
formatSpec = '%f';
mec = fscanf(mecFileId, formatSpec);
nRows = 3; % wavelength, HbO2, Hb
nCols = length(mec) / nRows;
mec = reshape(mec, [nRows, nCols]);
lambdas = mec(1, :);
mec_HbO2 = mec(2, :);
mec_Hb = mec(3, :);

lambda_fixed = lambdas(1);
assert(lambda_fixed == 532);

mecMat = zeros(2,2);
mecMat(1,1) = mec_HbO2(1);
mecMat(1,2) = mec_Hb(1);

nLambda = length(lambdas);
det_list = zeros(nLambda-1, 1);
condiNumber_list = zeros(nLambda-1, 1);

for iLambda = 2 : nLambda

    mecMat(2,1) = mec_HbO2(iLambda);
    mecMat(2,2) = mec_Hb(iLambda);
    
    currDet = det(mecMat);
    det_list(iLambda-1) = currDet;
    
    currCondiNumber = cond(mecMat);
    condiNumber_list(iLambda-1) = currCondiNumber;
    
end

figure; 
subplot(211); plot(lambdas(2:end), condiNumber_list);
xlabel('Wavelengths'); ylabel('Condition number');
xlim([lambdas(2), lambdas(end)]);
title('Condition number plot');
subplot(212); plot(lambdas(2:end), det_list);
hold on; plot(lambdas(2:end), zeros(nLambda-1,1), 'r--'); hold off;
xlabel('Wavelengths'); ylabel('Determinant');
xlim([lambdas(2), lambdas(end)]);
title('Determinant plot');




