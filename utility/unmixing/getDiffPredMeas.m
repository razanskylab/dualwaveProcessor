function diffPredMeas = getDiffPredMeas(sO2, fwdModel, measSig, ...
                            zVoxel_mm, zVoxelSize_mm, vesselLoc_idx, ...
                            mecMat, flagDiffMethod)
% calculate the "difference" between model prediction and measurement as a
% function of sO2 level
% the exact "difference" definition depends on the input "flagDiffMethod"

nWavelength = size(measSig, 2);
vesselStartIdx = vesselLoc_idx(1);
vesselEndIdx = vesselLoc_idx(2);
nVesselIdx = length(vesselStartIdx:vesselEndIdx);

% sO2 -> miu_a
miu_a_vec = zeros(length(zVoxel_mm), nWavelength); % [nz,nWavelength]
% mecMatNrmConst = max(mecMat(:));
% mecMatNrm = mecMat ./ mecMatNrmConst;
miu_a = mecMat * [sO2; (1-sO2)]; % [nWavelength,1]
% [~, miuaStartIdx] = min(abs(zVoxel_mm - zVessel_start));
% [~, miuaEndIdx] = min(abs(zVoxel_mm - zVessel_end));
miu_a_vec(vesselStartIdx:vesselEndIdx, :) = repmat(miu_a', nVesselIdx, 1);
% miu_a -> phi
phi_vec = exp(-(cumsum(miu_a_vec).*(zVoxelSize_mm*1e-1))); % [nz,nWavelength]
% miu_a_nrm = miu_a ./ mecMatNrmConst; % [nz,nWavelength]
% phi -> p0
p0_vec = miu_a_vec .* phi_vec; % [nz,nWavelength]

if nWavelength == 1
    if strcmp(flagDiffMethod, 'L2NormSquared')
        Ax = fwdModel * p0_vec;
        diffPredMeas = norm(Ax-measSig)^2;
    else
        error('Invalid difference method for single wavelength case!');
    end
elseif nWavelength == 2
        Ax_lambda1 = squeeze(fwdModel(:,:,1)) * p0_vec(:,1); % [nt,1], 532
        Ax_lambda2 = squeeze(fwdModel(:,:,2)) * p0_vec(:,2); % [nt,1], 558
        Ax = [Ax_lambda1, Ax_lambda2]; % [nt, 2]
        ampPred = max(Ax, [], 1); % [1,2]
        % compute normalized spectral slope
        dcOffsetPred = mean(ampPred);
        ampPredNrm = ampPred ./ dcOffsetPred;
        predSpectralSlope = diff(ampPredNrm); % ampPredNrm(2)-ampPredNrm(1)
        ampMeas = max(measSig, [], 1); % [1,2]
        dcOffsetMeas = mean(ampMeas);
        ampMeasNrm = ampMeas ./ dcOffsetMeas;
        measSpectralSlope = diff(ampMeasNrm); % ampMeasNrm(2)-ampMeasNrm(1)
    if strcmp(flagDiffMethod, 'L2NormSquared')
        diffPredMeas = norm(Ax-measSig)^2;
    elseif strcmp(flagDiffMethod, 'spectralSlope')
        diffPredMeas = (predSpectralSlope - measSpectralSlope)^2;
    elseif strcmp(flagDiffMethod, 'amplitudeRatio')
        diffPredMeas = ((ampPred(2)/ampPred(1))-(ampMeas(2)/ampMeas(1)))^2;
    else
        error('Invalid difference method in dual-wavelength case!');
    end
else
    error('Number of wavelengths can ONLY be 1 or 2!');
end

end
