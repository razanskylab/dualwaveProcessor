function reduScatterCoeff = get_reduced_scattering_coeff_skin(lambda)
%Calculates reduced scattering coefficient in the skin (epidermis & dermis)
%with the empirical formula described in the following paper:
% Effects of wavelength-dependent fluence attenuation on the noninvasive
% photoacoustic imaging of hemoglobin oxygen saturation in subcutaneous 
% vasculature in vivo, Inverse Problems 23 (2007) S113-S122.

% INPUT: lambda, wavelength in [nm]
% OUTPUT: reduScatterCoeff, reduced scattering coefficient in [cm-1]

reduScatterCoeff = 1.1*((10^12)*(1/(lambda^4))) + 73.7*(1/(lambda^0.22));

end