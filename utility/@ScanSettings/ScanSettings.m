%% File: ScanSettings.m @ ScanSettings
%  Date: 30.06.2021
%  Author: Urs Hofmann, Weiye Li - RazanskyLab

classdef ScanSettings < handle

	properties
		scanName(1,:) char = 'new scan';
		dataBase(1,:) char = 'C:\data_temp\'; % parent directory to store scan data
		outFolder(1,:) char = 'C:\data_temp\today_date'; % directory pointing to today's folder
		outFileName(1,:) char = '001_test'; % output file name without suffix (can be .mat or .jpg)
		outPathFull(1,:) char = 'C:\data_temp\today_date\001_test.mat'; % full directory to store scan data (up to .mat)

		flagSaveRawData(1,1) {mustBeNumericOrLogical} = 1;
		flagLivePreview(1,1) {mustBeNumericOrLogical} = 0;

		% field-of-view
		% x direction is fast/VoiceCoil stage
		% y direction is slow/linear Thorlabs stage
		width(1,2) {mustBeNumeric, mustBeFinite} = [10, 10]; % [mm], [x,y]
		dr(1,2) {mustBeNumeric, mustBeFinite} = [20, 20] .* 1e-3; % [mm], [x,y]
		center(1,2) {mustBeNumeric, mustBeFinite} = [6, 25]; % [mm], [x,y]

		% stage
		velY(1,1) {mustBeNumeric, mustBeFinite} = 20; % [mm/s], velocity of Y(slow/linear) stage
		maxVelX(1,1) {mustBeNumeric, mustBeFinite} = 200; % [mm/s], max. velocity of X stage
		
		% DAQ settings
		nSamples(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 2048;
		sensitivityUs(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 1000; % [mV]
		sensitivityPd(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 5000; % [mV]
		samplingFreq(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 250e6; % [Hz]
		nAverages(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 1;
		% cropping indices along A-scan 
		pdCrop(1,2) {mustBeInteger, mustBeNumeric, mustBeFinite} = [1, 1000];
		usCrop(1,2) {mustBeInteger, mustBeNumeric, mustBeFinite} = [1, 1024];
		pdNoiseWindow(1,2) {mustBeInteger, mustBeNumeric, mustBeFinite} = [1, 50];

		% laser parameters
		wavelengths(1,:) {mustBeInteger, mustBeNumeric, mustBeFinite} = [532, 558]; % [nm]
		maxPrf(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 10e3; % [Hz]

		% used for calculating depth vector during live preview
		pulserDelay = 0.6; % [us], delay time of US pulser
		pdInitPeak_532(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 150; % [index]
		pdInitPeak_dye(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 220; % [index]

		nColors(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 40;
		nDepthLabels(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite} = 15;
		transparency(1,1) {mustBeNumeric, mustBeFinite} = 1; % can be e.g. [0.5, 2]
		minAmp(1,1) {mustBeNumeric, mustBeFinite} = 0.2; % noise floor in normalized MIP considered during live preview

		SOS(1,1) {mustBeNumeric, mustBeFinite} = 1500; % [m/s], speed-of-sound
		fd(1,1) single = 7e-3; % focal distance of transducer
		bScanRate(1,1) {mustBeNumeric, mustBeFinite}; % automatically calculated, do not set me

		% missedTriggerCount(1,1) double = 0; % keep track of missed trigger in case it happens, for later B-scan shift correction

	end

	properties(Dependent)
		nWavelengths(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		nUs(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		nPd(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		xVec(1,:) {mustBeNumeric, mustBeFinite};
		yVec(1,:) {mustBeNumeric, mustBeFinite};
		tVec(1,:) {mustBeNumeric, mustBeFinite};
		zVec(2,:) {mustBeNumeric, mustBeFinite};
		% contains dual-wave z vectors, first row is pulse-echo / onda532, second row is Dye
		originPoint(1,2) {mustBeNumeric, mustBeFinite};
		endPoint(1,2) {mustBeNumeric, mustBeFinite};
		nX(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		nY(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		nBScans(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		nALines(1,1) {mustBeInteger, mustBeNumeric, mustBeFinite};
		dt(1,1) {mustBeNumeric, mustBeFinite};
		mode(1,:) char;
	end

	methods

		% calculate dependent properties with get functions

		function mode = get.mode(ss)
			if length(ss.wavelengths) == 1
				if ss.wavelengths == 0
					mode = 'us';
				elseif ss.wavelengths == 532
					mode = 'onda532';
				elseif ss.wavelengths > 532
					mode = 'dye';
				else
					error('Unknown mode at single wavelength!');
				end
			elseif length(ss.wavelengths) == 2
				mode = 'dualwave';
			else
				error('Length of wavelengths MUST be either 1 or 2!');
			end
		end

		function nBScans = get.nBScans(ss)
			nBScans = ss.width(2) / ss.dr(2);
			if mod(nBScans, 2)
				error('[VoiceScan] # of B-scans must be even!');
			end
		end

		function nALines = get.nALines(ss)
			nALines = ss.nX;
		end

		function nWavelengths = get.nWavelengths(ss)
			nWavelengths = length(ss.wavelengths);
		end

		function dt = get.dt(ss)
			dt = 1 / ss.samplingFreq;
		end

		function nPd = get.nPd(ss)
			nPd = ss.pdCrop(2) - ss.pdCrop(1) + 1;
		end

		function nUs = get.nUs(ss)
			nUs = ss.usCrop(2) - ss.usCrop(1) + 1;
		end

		function originPoint = get.originPoint(ss)
			originPoint = ss.center - ss.width./2;
		end

		function endPoint = get.endPoint(ss)
			endPoint = ss.center + ss.width./2;
		end

		function xVec = get.xVec(ss)
			% compute x vectors in [mm]
			xVec = ss.originPoint(1) : ss.dr(1) : ss.endPoint(1);
		end

		function yVec = get.yVec(ss)
			% compute y vectors in [mm]
			yVec = ss.originPoint(2) : ss.dr(2) : ss.endPoint(2);
		end

		function tVec = get.tVec(ss)
			% computes time vector in [s] without cropping
			tVec = (0:(ss.nSamples-1)) .* ss.dt;
		end

		function zVec = get.zVec(ss)
			% compute z (depth) vectors in [mm] with cropping
			usCrop = ss.usCrop;
			usCropIndices = (usCrop(1):usCrop(2));
			nUs = length(usCropIndices);
			zVec = zeros(2, nUs);
			if strcmp(ss.mode, 'us')
				% US pulse-echo
				usInitPeak = round(ss.pulserDelay * 1e-6 * ss.samplingFreq); % 150 if sampling @ 250 MHz
				t_samples = usCropIndices - usInitPeak + 1;
				t = t_samples .* ss.dt;
				zVec(1,:) = t .* ss.SOS ./ 2 .* 1e3 - ss.fd * 1e3; % [1, nUs]
			elseif strcmp(ss.mode, 'onda532')
				% Onda532 only
				t_samples = usCropIndices - ss.pdInitPeak_532 + 1;
				t = t_samples .* ss.dt;
				zVec(1,:) = t .* ss.SOS .* 1e3  - ss.fd * 1e3; % [1, nUs]
			elseif strcmp(ss.mode, 'dye')
				% Dye laser only
				t_samples = usCropIndices - ss.pdInitPeak_dye + 1;
				t = t_samples .* ss.dt;
				zVec(2,:) = t .* ss.SOS .* 1e3  - ss.fd * 1e3; % [1, nUs]
			elseif strcmp(ss.mode, 'dualwave')
				% dual-wavelength
				t_samples = usCropIndices - ss.pdInitPeak_532 + 1;
				t = t_samples .* ss.dt;
				zVec(1,:) = t .* ss.SOS .* 1e3  - ss.fd * 1e3;
				t_samples = usCropIndices - ss.pdInitPeak_dye + 1;
				t = t_samples .* ss.dt;
				zVec(2,:) = t .* ss.SOS .* 1e3  - ss.fd * 1e3; % [2, nUs]
			else
				error('[VoiceScan] Unknown scan mode when calculating depth vector!');
			end

		end

		function nX = get.nX(ss)
			nX = length(ss.xVec);
		end

		function nY = get.nY(ss)
			nY = length(ss.yVec);
		end

		% calculate public properties with set functions
		function set.width(ss, width)
			ss.width = width;
		end

		function set.dr(ss, dr)
			ss.dr = dr;
		end

		function set.center(ss, center)
			ss.center = center;
		end

		function set.nAverages(ss, nAverages)
			if (nAverages < 1)
				error("Number of averages must be at least 1");
			end
			ss.nAverages = nAverages;
		end

		function set.wavelengths(ss, wavelengths)

			% check if number of wavelengths is appropriate
			if length(wavelengths) > 2
				error("Cannot use more then two wavelengths");

			elseif length(wavelengths) == 2
				if any(wavelengths < 532)
					error("Cannot fire below 532 nm!");
				end
				if any(wavelengths > 720)	
					error("Cannot fire above 720 nm!");
				end
				if (wavelengths(1) == wavelengths(2))
					error("Wavelengths cannot be the same");
				end
				if ~(any(wavelengths == 532))
					error("In dual wavelength mode, one needs to be 532");
				end
				ss.wavelengths = wavelengths;

			elseif length(wavelengths) == 1
				if (wavelengths(1) ~= 0) && (wavelengths(1) < 532)
					error("Cannot fire below 532 nm!");
				end
				if (wavelengths(1) ~=0) && (wavelengths(1) >720)
					error("Cannot fire above 720 nm!");
				end
				ss.wavelengths = wavelengths;

			else
				error('MUST specify wavelengths!');

			end

		end

		Validate_Settings(ss);

	end

end