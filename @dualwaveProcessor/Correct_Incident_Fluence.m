% File: Correct_Incident_Fluence.m @ dualwaveProcessor
% Date: July 28th, 2022
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

% Change log on Feb. 19, 2023: Disabled scaling the whole US signal volume with
% ppeCorrectionMap. Moved scaling to the unmixing function.

% Major update date: August 23, 2023
% Log: changed function name from Correct_PPE_Fluctuation to Correct_Incident_Fluence
% Current strategy: first normalize within wavelength to the mean/median PPE level, then normalize between wavelength

function Correct_Incident_Fluence(dwp)

    fprintf('[dualwaveProcessor] Correcting incident fluence variation within and between wavelengths...');

    if strcmp(dwp.scanSett.mode, 'dualwave')

      % covert PD raw signal into PPE estimation in arbitrary unit
      dwp.ppeMap_532 = pd_shot_to_energy_au(dwp.pdVol_532); % [nx, ny]
      dwp.ppeMap_dye = pd_shot_to_energy_au(dwp.pdVol_dye); % [nx, ny]
      ppeMap_532 = dwp.ppeMap_532;
      ppeMap_dye = dwp.ppeMap_dye;

      % two-step correction to bring PD reading at dye wavelength to the actual PPE at the output
      % ONLY care about estimating the true relative PPE level between 532 and 558
      % step 1: in-fiber PD responds to dye wavelengths more than 532nm
      if dwp.flagPDResponseScaling
        ppeMap_dye = ppeMap_dye ./ dwp.pdResonsivityRatio;
      else
        warning('[dualwaveProcessor] PD response correction DISABLED!');
      end
      % step 2: splitter fiber passes dye wavelengths more than 532nm
      if dwp.flagFiberSplitScaling
        ppeMap_dye = ppeMap_dye ./ dwp.fiberSplitRatio;
      else
        warning('[dualwaveProcessor] fiber splitting ratio correction DISABLED!');
      end

      if dwp.flagPPEMask
        nStd = 5; % assuming Gaussian distribution of PPEs, how many std to consider as good ones
        mean532 = mean(ppeMap_532(:), 'omitnan');
        std532 = std(ppeMap_532(:), 'omitnan');
        meanDye = mean(ppeMap_dye(:), 'omitnan');
        stdDye = std(ppeMap_dye(:), 'omitnan');
        ppeLB532 = mean532 - nStd*std532;
        ppeUB532 = mean532 + nStd*std532;
        ppeLBDye = meanDye - nStd*stdDye;
        ppeUBDye = meanDye + nStd*stdDye;
        ppeMask_532 = (ppeMap_532>ppeLB532) & (ppeMap_532<ppeUB532);
        ppeMask_dye = (ppeMap_dye>ppeLBDye) & (ppeMap_dye<ppeUBDye);
        dwp.ppeMask = single(ppeMask_532 & ppeMask_dye);
        % already mask out outlier PPEs
        ppeMap_532 = ppeMap_532 .* dwp.ppeMask;
        ppeMap_dye = ppeMap_dye .* dwp.ppeMask;
        ppeMap_532(ppeMap_532==0) = nan;
        ppeMap_dye(ppeMap_dye==0) = nan;
      else
        warning('[dualwaveProcessor] outlier PPE thresholding DISABLED!');
      end

      % visualize PPE histograms
      if dwp.flagShowPpeHist
        allPpes = [ppeMap_532(:); ppeMap_dye(:)];
        maxPpe = max(allPpes);
        minPpe = min(allPpes);
        binStepSize = 100;
        binEdges = minPpe:binStepSize:maxPpe;
        fHist = figure('Name', 'Histogram of PPEs', 'Color', 'w');
        % fHist.Units = 'normalized';
        % fHist.OuterPosition = [0 0 1 1];
        histogram(ppeMap_532(:), binEdges);
        hold on; histogram(ppeMap_dye(:), binEdges);
        hold off; grid on;
        legend('532', '558'); title('PPEs [a.u.]');
      end

      % compute within-wavelength & between-wavelength scaling factors
      if strcmp(dwp.refPpeMethod, 'mean')
        % within wavelength
        dwp.ppeCorrectionMap_532 = mean(ppeMap_532(:),'omitnan') ./ ppeMap_532;
        dwp.ppeCorrectionMap_dye = mean(ppeMap_dye(:),'omitnan') ./ ppeMap_dye;
        % between wavelength
        dwp.scaleFctIncident = mean(ppeMap_532(:),'omitnan') / mean(ppeMap_dye(:),'omitnan'); 
      elseif strcmp(dwp.refPpeMethod, 'median')
        % within wavelength
        dwp.ppeCorrectionMap_532 = median(ppeMap_532(:),'omitnan') ./ ppeMap_532;
        dwp.ppeCorrectionMap_dye = median(ppeMap_dye(:),'omitnan') ./ ppeMap_dye;
        % between wavelength
        dwp.scaleFctIncident = median(ppeMap_532(:),'omitnan') / median(ppeMap_dye(:),'omitnan');
      elseif strcmp(dwp.refPpeMethod, 'mode')
        % within wavelength
        dwp.ppeCorrectionMap_532 = mode(ppeMap_532(:)) ./ ppeMap_532;
        dwp.ppeCorrectionMap_dye = mode(ppeMap_dye(:)) ./ ppeMap_dye;
        % between wavelength
        dwp.scaleFctIncident = mode(ppeMap_532(:)) / mode(ppeMap_dye(:));
      else
        error('[dualwaveProcessor] invalid method to calculate reference PPE level!');
      end

    else
      error('[dualwaveProcessor] PPE correction ONLY implemented for dualwave mode!');
    end

    fprintf('done!\n');

end