% File: correct_laser_jitter.m
% Weiye Li - Razansky Lab - University and ETH Zurich
% E-mail: weiye.li@uzh.ch
% Date: August 18, 2022

function [pdSig, usSig, pdMaxIdx] = correct_laser_jitter(pdSig, usSig)
%Corrects for laser jitter in PD and US signals by shifting each
%A-line w.r.t. a reference PD signal temporal position (taken as the
%median of the peak location of all recorded PD signals)

% INPUT: uncorrected PD and US signals, can be 2D or 3D arrays, first dimension is time (t) 
% OUTPUT: corrected PD and US signals, same dimension as the inputs, pdMaxIdx is the reference temporal position of PD signals

% check dimensions of pd and us signals
nDim = numel(size(pdSig));
if nDim ~= numel(size(usSig))
    error('US and PD signals MUST have the same number of dimensions! Can be 2 or 3 dimensional arrays.');
end

% reshape into 2D matrices if inputs are 3D volumes
if nDim == 3
    [ntPd, nxPd, nyPd] = size(pdSig);
    [ntUs, nxUs, nyUs] = size(usSig);
    assert(nxPd == nxUs);
    assert(nyPd == nyUs);
    pdSig = reshape(pdSig, [ntPd, nxPd*nyPd]);
    usSig = reshape(usSig, [ntUs, nxUs*nyUs]);
end

% get median PD max index and calculate how much to shift for each shot
[~, pdMaxIndices] = max(pdSig, [], 1);
pdMaxIdx = median(pdMaxIndices);
pdShiftIndices = pdMaxIdx - pdMaxIndices;

% loop over each scan position
for ixy = 1 : size(pdSig, 2)
    currShiftIdx = pdShiftIndices(ixy);
    if (currShiftIdx == 0)
        continue;
    else
        pdSig(:,ixy) = circshift(pdSig(:,ixy), currShiftIdx);
        usSig(:,ixy) = circshift(usSig(:,ixy), currShiftIdx);
    end
end

% reshape back to 3D volumes if needed
if nDim == 3
    pdSig = reshape(pdSig, [ntPd, nxPd, nyPd]);
    usSig = reshape(usSig, [ntUs, nxUs, nyUs]);
end

end