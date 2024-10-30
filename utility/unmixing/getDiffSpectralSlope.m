function diffSpectralSlope = getDiffSpectralSlope(miu_a, fwdModel, measSpectralSlope, voxelSize, mecMatNrmConst)
% calculate the "normalized" spectral slope difference between model
% prediction and measured signal, normalization here means to get rid of
% the "DC offset" of the amplitude spectrum between prediction and
% measurement

% independent of # of wavelengths
phi = exp(-(cumsum(miu_a).*(voxelSize*1e2))); % [nz,2], voxelSize is in [m]
miu_a_nrm = miu_a ./ mecMatNrmConst; % [nz,2]
p0 = miu_a_nrm .* phi; % [nz,2]

nWavelength = size(miu_a, 2);
if nWavelength == 1
    error('Number of wavelengths MUST be at least 2 to calculate spectral slope!');
elseif nWavelength == 2
    % compute model prediction
    AxLambda1 = squeeze(fwdModel(:,:,1)) * p0(:,1); % [nt,1], 532
    AxLambda2 = squeeze(fwdModel(:,:,2)) * p0(:,2); % [nt,1], 558
    % compute normalized spectral slope
    ampLambda1 = max(AxLambda1);
    ampLambda2 = max(AxLambda2);
    dcOffset = mean([ampLambda1, ampLambda2]);
    ampSpectraNrm = [ampLambda1, ampLambda2] ./ dcOffset;
    predSpectralSlope = ampSpectraNrm(2) - ampSpectraNrm(1);
    diffSpectralSlope = (predSpectralSlope - measSpectralSlope)^2;
else
    error('Number of wavelengths can NOT be larger than 2!');
end

end
