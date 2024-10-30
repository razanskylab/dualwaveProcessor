function Validate_Settings(ss)

	% check if we do not crash into bouandries of stages
	if any((ss.center - (ss.width / 2)) < 0.5)
		error("Out of stage movement range!");
	end

	upLim = [12, 50]; % max. movement range of stage
	if any((ss.center + (ss.width / 2)) > upLim)
		error("Out of stage movement range!");
	end

	if (ss.usCrop(2) < ss.usCrop(1))
		error("Invalid cropping of us channel!");
	end


	if (ss.pdCrop(2) < ss.pdCrop(1))
		error("Invalid cropping of pd channel!");
	end

	if any([ss.usCrop(2), ss.pdCrop(2)] > ss.nSamples)
		error("Cropping indices exceed total number of samples!");
	end

end