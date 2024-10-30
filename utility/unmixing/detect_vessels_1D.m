function [nVessel, vesselIndices] = detect_vessels_1D(rawAline, vesselSigThreshold, maxNVessels)
%Perform vessel detection along depth dimention
% ------- Strategy -------
% identify # of detected max. amplitude peaks that are at least
% axial resolution apart          
% i) if only one max. amplitude peak, meaning only one vessel, we then
%    estimate the vessel size by the distance between the positive peak
%    and the immediate-after negative peak;
% ii) if two max. amplitude peaks detected, meaning two vessels,
%     (we assume only two layers of vessels are excitable),
%     we estimate their respective size as above.

% INPUT:
% - rawAline: raw 1D A-line signal, [nt,1]
% - vesselSigThreshold: the signal amplitude threshold above which is
% considered legit vessel signal
% - maxVesselSize_idx: max. possible vessel size in unit of indices
% - axialResolu_idx: axial resolution in unit of indices

% OUTPUT:
% - nVessel: # of vessels detected, can be either 0, 1, or 2
% - vesselIndices: the indices of detected vessel peaks
% - vesselSizes_idx: the estimated sizes of detected vessels in [indices]

nVessel = 0;
vesselIndices = zeros(1,maxNVessels);
% vesselSizes_idx = zeros(1,maxNVessels);
for iVessel = 1 : maxNVessels
    [maxVal,maxIdx] = max(rawAline);
    if maxVal < vesselSigThreshold
        % NO legit vessel signal detected
        break;
    else
        % detected a legit vessel signal
        vesselIndices(iVessel) = maxIdx;
%         searchIdxEnd = min([maxIdx+maxVesselSize_idx, nt]);
%         searchSegment = rawAline(maxIdx+1 : searchIdxEnd);
%         [~,vesselEndIdx] = min(searchSegment);
%         vesselSizes_idx(iVessel) = vesselEndIdx;
        % after extracting info. about the current vessel, remove it to
        % prepare for the detection of the 2nd vessel
%         i1 = max([(maxIdx-axialResolu_idx), 1]);
%         i2 = min([(maxIdx+vesselEndIdx+axialResolu_idx), nt]);
%         rawAline(i1:i2) = vesselSigThreshold/2;
    end
    
    nVessel = nVessel + 1;
        
end

end

