% File: dualwaveProcessor.m @ dualwaveProcessor
% Date: July 28th, 2022
% Author: Weiye Li, Razansky Lab, University and ETH Zurich
% E-mail: weiye.li@uzh.ch

% ------- Major update on August 23, 2023 -------
% Log: implemented full processing pipeline from raw data to 3D sO2 image
% Main function blocks are:
% 1) Basic pre-processing - like averaging, correcting scan pattern, coverting data type, depth cropping
% 2) Incident fluence correction - both in-wavelength and between-wavelength normalization
% 3) Skin/skull surface estimation with US
% 4) Local fluence correction - wavelength-dependent fluence variation after light-tissue interaction
% 5) Vessel reconstruction along depth - estimate vessel size along A-line (unmix, then reconstruct)
% 6) Render 3D sO2 image from reconstructed vessels and estimated sO2 - assign sO2 value to vessel segment along depth, visualize in 3D

classdef dualwaveProcessor < handle

  properties

    data_dir(1,:) char = 'Please define data folder!';
    data_name(1,:) char = 'Please define data name!';

    psfFile_dir(1,:) char = 'C:\data_temp\2023_08_17\sphereModel_dw_processed.mat';

    molarExtinctCoeff_dir(1,:) char = 'C:\Code\dualwave-foam-dev\09_Utility\unmixing';
    molarExtinctCoeff_fileName(1,:) char = 'Molar_Extinction_Coefficients_Hemoglobin_Water_532_to_610.txt';

    scanSett(1,1) ScanSettings;

    RawDataPd(:,:,:,:,:) int16; % [t, lambda, averages, x, y]
    RawDataUs(:,:,:,:,:) int16; % [t, lambda, averages, x, y]

    usVol_532(:,:,:) single; % [t, x, y]
    usVol_dye(:,:,:) single; % [t, x, y]
    usVol_pe(:,:,:) single; % [t, x, y], pulse-echo volume

    usVolTissue_532(:,:,:) single; % [nPixelTissue, nx, ny], extracted tissue data
    usVolTissue_dye(:,:,:) single; % [nPixelTissue, nx, ny], extracted tissue data
    sO2VolTissue(:,:,:) single; % [nPixelTissue, nx, ny], extracted tissue data

    pdVol_532(:,:,:) single; % [t, x, y]
    pdVol_dye(:,:,:) single; % [t, x, y]

    ppeMap_532(:,:) single; % [x, y]
    ppeMap_dye(:,:) single; % [x, y]

    ppeCorrectionMap_532(:,:) single; % [nx, ny]
    ppeCorrectionMap_dye(:,:) single; % [nx, ny]

    psf532(:,:) single; % [nt, nz]
    psfDye(:,:) single; % [nt, nz]
    psfVoxel_mm(1,:) single; % [1, nz]
    psfVec_mm_532(1,:) single; % [1, nt]
    psfVec_mm_dye(1,:) single; % [1, nt]

    snrMask(:,:) single; % [x, y], used together with "flagSNRMask", only consider pixels with high-enough intensity
    ppeMask(:,:) single; % [x, y], get rid of outlier PPEs

    sO2_map(:,:) single; % [x, y]
    vesselDepth_map(:,:) single; % [x, y]
    vesselStructure_map(:,:) single; % [x, y]
    sO2_vol(:,:,:) single; % [z, x, y]
    attenuation_map(:,:) single; % [x, y]

    pdResonsivityRatio(1,1) single = 1.1; % Dye / Onda, based on Thorlabs plot
    fiberSplitRatio(1,1) single = 1.09; % Dye / Onda, based on Thorlabs plot
    % Calibrated product of the above two ratios is ~1.2

    refPpeMethod(1,:) char = 'median';

    daqResolution(1,1) single = 16;

    unmix2DMethod(1,:) char = 'peak2peak'; % 'max', 'peak2peak'

    assignToTissueLayerInfo(1,:) char = 'structure'; % 'structure', 'sO2'

    sO2_list(1,:) single = 0 : 0.1 : 1; % list of physiologically reasonable sO2 levels

    depthStart(1,1) single = 6e-3; % [m], depth coverage closest to transducer
    depthEnd(1,1) single = 9e-3; % [m], depth coverage farest away from transducer
    % these are set to be +/- 1.5 mm around the focal distance (~7.5 mm)
    depthStartRecon(1,1) single = 7e-3; % [m], starting depth for reconstruction z voxel grid
    depthEndRecon(1,1) single = 8e-3; % [m], end depth for reconstruction z voxel grid
    zVoxelSizeRecon(1,1) single = 5e-6; % [m], voxel size of reconstruction z voxel grid

    flagCorrectPpeFluctuation(1,1) logical = true;
    flagPDResponseScaling(1,1) logical = true;
    flagFiberSplitScaling(1,1) logical = true;
    flagShowPpeHist(1,1) logical = true;
    flagSNRMask(1,1) logical = true; % whether to apply SNR thresholding to show unmixed image
    flagPPEMask(1,1) logical = true; % whether to apply PPE thresholding to remove outlier PPEs (out of +/- 2*sigma range)
    
    % Normally activated in unmixing mode to have more stable PPEs at the center of FOV
    flagCropFOVBoundary(1,1) logical = false; % whether to crop out starting and end movement range of fast stage
    cropXBoundary_mm = 0.5;
    cropYBoundary_mm = 0.5;

    flagZoominScan(1,1) logical = false; % whether this is a zoom-in scan
    fullFOVCenter(1,2) double = [6, 25]; % [x,y] center position of scanning stages when acquiring full FOV images
    zoominFOVCenter(1,2) double = [6, 25]; % [x,y] center position of scanning stages when acquiring zoom-in images

    amplitudeLowerThresholdPct(1,1) single = 0.05;
    amplitudeUpperThresholdPct(1,1) single = 0.75;

    stepSizeWound_mm_us(1,1) single = 10e-3;
    stepSizeWound_mm_oa(1,1) single = 5e-3;
    stepSize_mm_zoomin(1,1) single = 2e-3;

    tissueCropDepths_mm_us(:,:) single;
    tissueCropDepths_mm_oa(:,:) single;

    aboveTissue_mm(1,1) single = 0.018;
    belowTissue_mm(1,1) single = 0.72;

    epidermisThickness_mm(1,1) single = 0.018;
    dermisThickness_mm(1,1) single = 0.198;
    junctionThickness_mm(1,1) single = 0.1;

    vesselSigLB(1,1) single = 10; % [mV], lower bound of vessel signal intensity
    maxNVessels(1,1) single = 1; % assume maximum one vessel per A-line

    scaleFctIncident(1,1) single = 1; % between wavelength incident PPE level scaling factor, to be multiplied with Dye dataset

    Reflectance532(1,1) single = 0.972;
    ReflectanceDye(1,1) single = 0.868;
    % These values are based on visual estimates of Fig. 3 data in the following paper: 
    % Gender variations in the optical properties of skin in murine animal models,
    % Journal of Biomedical Optics, Vol. 16, Issue 1, 011008 (January 2011)
    % Assumed linear model between these two points: (525, 1) and (575, 0.8)
    % then the line equation is: R = -0.004*lambda + 3.1

    Cm(1,1) single = 3; % [%]
    % assumed volume fraction of melanosome in the epidermis, varies from 1.3% to 43% depending on skin color
    % Refer to the following paper:
    % Effects of wavelength-dependent fluence attenuation on the noninvasive
    % photoacoustic imaging of hemoglobin oxygen saturation in subcutaneous 
    % vasculature in vivo, Inverse Problems 23 (2007) S113-S122.

    % flagSO2Threshold(1,1) logical = false; % whether to threshold sO2 map to [0,1]
    % flagSpectralSlopeMask(1,1) logical = false;

  end

  properties(Dependent)

    % cropping indices along depth depends on desired depth coverage and focal distance
    zCropIdx_us(2,2) single = ones(2,2); % first row is for pulse-echo / 532, second row is for Dye
    zVec_mm(2,:) single; % [nWavelength, nt], first index is 532, 2nd is 558
    zVoxel_mm(1,:) single; % [1, nt], Z voxel grid (on which to reconstruct) is the same for both wavelengths

    % molar extinction coefficient matrix
    mecMat(2,2) single; % [HbO2(532), Hb(532);
                        %  HbO2(Dye), Hb(Dye)].

    % reference blood absorption spectra, i.e. blood absorption coefficient as a function of wavelength at different sO2 levels
    miu_a_mat(2,:) single;

    spectralSlope_list(1,:) single;

    samplingDist_mm(1,1) single;

    scaleFctReflectance(1,1) single;

  end

  methods

    % MAIN functionalities
    Load_Raw_Data(dwp); % load 5D raw data, already crop in depth based on "zCropIdx_us"
    Preproc_Basic(dwp); % correct jitter, correct raster scan pattern, convert count to voltage, average
    Correct_Incident_Fluence(dwp); % normalize within- and between-wavelength fluence variations
    Load_Tissue_Surface(dwp); % load estimated skin/skull surface from US dataset and interpolate onto OA scanning grid
    Unmix_2D(dwp); % unmix only using max. amplitudes
    Unmix_3D(dwp); % unmix while considering wavelength-dependent fluence attenuation along depth
    Render_sO2_3D(dwp); % render 3D structure and sO2 image, enable visualizations including layered top view, side view and volumetric view
    Extract_Tissue_Data(dwp); % after all the processing, extract only useful data from the tissue, basically depth cropping based on skin surface
    [sO2_map_dermis, sO2_map_hypodermis] = Assign_To_Tissue_Layer(dwp);
    [structure_map_junction, sO2_map_junction] = Visualize_Junction_Vasculature(dwp);

    % BACKUP functionalities
    Load_PSF(dwp); % load measured and processed sphere model for vessel reconstruction
    Reconstruct_Vessel_1D(dwp); % one-dimensional vessel detection and reconstruction along depth
    
    % TODO:
    % Estimate_Tissue_Surface(dwp); % interactive function to delineate skin/skull surface based on US dataset
    % % FOR NOW get the tissue surface data in a separate script and just load it into dualwaveProcessor
  
  end

  methods % getters
  
    function zCropIdx_us = get.zCropIdx_us(dwp)
      z_mm = dwp.scanSett.zVec + dwp.scanSett.fd*1e3; % [2, nZ], each row corresponds to one wavelength
      mode = dwp.scanSett.mode; % 'dualwave', 'onda532', 'dye', 'us'

      if strcmp(mode, 'dualwave')
        [~, idxStart_532] = min(abs(z_mm(1,:) - dwp.depthStart*1e3));
        [~, idxEnd_532] = min(abs(z_mm(1,:) - dwp.depthEnd*1e3));
        [~, idxStart_dye] = min(abs(z_mm(2,:) - dwp.depthStart*1e3));
        [~, idxEnd_dye] = min(abs(z_mm(2,:) - dwp.depthEnd*1e3));
        zCropIdx_us(1,:) = single([idxStart_532, idxEnd_532]);
        zCropIdx_us(2,:) = single([idxStart_dye, idxEnd_dye]);

      elseif strcmp(mode, 'onda532')
        [~, idxStart] = min(abs(z_mm(1,:) - dwp.depthStart*1e3));
        [~, idxEnd] = min(abs(z_mm(1,:) - dwp.depthEnd*1e3));
        zCropIdx_us(1,:) = single([idxStart, idxEnd]);

      elseif strcmp(mode, 'dye')
        [~, idxStart] = min(abs(z_mm(2,:) - dwp.depthStart*1e3));
        [~, idxEnd] = min(abs(z_mm(2,:) - dwp.depthEnd*1e3));
        zCropIdx_us(2,:) = single([idxStart, idxEnd]);

      elseif strcmp(mode, 'us')
        [~, idxStart] = min(abs(z_mm(1,:) - dwp.depthStart*1e3));
        [~, idxEnd] = min(abs(z_mm(1,:) - dwp.depthEnd*1e3));
        zCropIdx_us(1,:) = single([idxStart, idxEnd]);

      else
        error('[dualwaveProcessor] invalid scan mode!');
      end

		end

    function mecMat = get.mecMat(dwp)
      % get molar extinction coefficients for current wavelength pair
      mecFileId = fopen(fullfile(dwp.molarExtinctCoeff_dir, dwp.molarExtinctCoeff_fileName), 'r');
      formatSpec = '%f';
      mec = fscanf(mecFileId, formatSpec);
      nRows = 3; % wavelength, HbO2, Hb
      nCols = length(mec) / nRows;
      mec = reshape(mec, [nRows, nCols]);
      lambdas = mec(1, :); mec_HbO2 = mec(2, :); mec_Hb = mec(3, :);
      mecMat = zeros(2, 2);
      idx_532 = find(lambdas == dwp.scanSett.wavelengths(1));
      idx_dye = find(lambdas == dwp.scanSett.wavelengths(2));
      mecMat(1,1) = mec_HbO2(idx_532);
      mecMat(1,2) = mec_Hb(idx_532);
      mecMat(2,1) = mec_HbO2(idx_dye);
      mecMat(2,2) = mec_Hb(idx_dye);
      mecMat = mecMat .* (150*2.303) ./ 64500; % now unit is [cm-1], see Prahl's page for details
    end

    function miu_a_mat = get.miu_a_mat(dwp)
      % compute reference blood absorption spectra
      mecMat = dwp.mecMat;
      sO2_list = dwp.sO2_list;
      miu_a_mat = zeros(2, length(sO2_list));
      for i = 1 : length(sO2_list)
        curr_sO2 = sO2_list(i);
        c_hbO2 = curr_sO2; c_hb = 1 - curr_sO2;
        miu_a_mat(:,i) = mecMat * [c_hbO2; c_hb];
      end
    end

    function spectralSlope_list = get.spectralSlope_list(dwp)
      miu_a_mat = dwp.miu_a_mat;
      wavelength_pair = dwp.scanSett.wavelengths;
      miu_a_mat = miu_a_mat - min(miu_a_mat, [], 1);
      miu_a_diff = miu_a_mat(2,:) - miu_a_mat(1,:);
      wavelength_diff = wavelength_pair(2) - wavelength_pair(1);
      spectralSlope_list = miu_a_diff ./ wavelength_diff;
    end

    function samplingDist_mm = get.samplingDist_mm(dwp)
      samplingDist_mm = (1/dwp.scanSett.samplingFreq) * dwp.scanSett.SOS * 1e3;
      samplingDist_mm = single(samplingDist_mm);
    end

    function scaleFctReflectance = get.scaleFctReflectance(dwp)
      scaleFctReflectance = (1/dwp.Reflectance532) / (1/dwp.ReflectanceDye);
    end

    function zVec_mm = get.zVec_mm(dwp)
      zVec_mm_532 = dwp.scanSett.zVec(1,:) + dwp.scanSett.fd*1e3;
      zVec_mm_532 = zVec_mm_532(dwp.zCropIdx_us(1,1):dwp.zCropIdx_us(1,2));
      zVec_mm_558 = dwp.scanSett.zVec(2,:) + dwp.scanSett.fd*1e3;
      zVec_mm_558 = zVec_mm_558(dwp.zCropIdx_us(2,1):dwp.zCropIdx_us(2,2));
      zVec_mm(1,:) = zVec_mm_532;
      zVec_mm(2,:) = zVec_mm_558;
    end

    function zVoxel_mm = get.zVoxel_mm(dwp)
      zVoxel_mm = (dwp.depthStartRecon : dwp.zVoxelSizeRecon : dwp.depthEndRecon) .* 1e3;
    end

  end

end