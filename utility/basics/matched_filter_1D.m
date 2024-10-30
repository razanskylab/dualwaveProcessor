function sig_mf = matched_filter_1D(Aline, psf1D, zVec_mm, zVoxel_mm, flagEnvelope)
%Performs matched filtering on an A-line with corresponding PSF

nt = length(zVec_mm);
nz = length(zVoxel_mm);
% make sure A-line and 1D PSF are column vectors with correction dimension
assert(length(Aline)==nt);
assert(length(psf1D)==nt);
Aline = reshape(Aline, [nt,1]);
psf1D = reshape(psf1D, [nt,1]);

% construct circulant matrix (convolution matrix) from PSF
Aconv = convmtx(psf1D, nz);

% NOW Aconv has row dimension of "full" convolution length, need to crop it
% to have length nt
[~,vesselDepthIdx] = max(Aline, [], 1);
vesselDepth_mm = zVec_mm(vesselDepthIdx);
[~,vesselVoxelIdx] = min(abs(zVoxel_mm - vesselDepth_mm));
% crop row dimension of conv. matrix to match nt
Aconv = Aconv(vesselVoxelIdx+1:vesselVoxelIdx+nt, :); % [nt,nz]

% matched filtering can be implemented by transpose multiplication
sig_mf = Aconv' * Aline;

% perform envelope extraction if enabled
if flagEnvelope
    [sig_mf, ~] = envelope(sig_mf);
end

end