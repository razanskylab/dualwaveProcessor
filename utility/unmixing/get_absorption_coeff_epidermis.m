function absorpCoeffEpi = get_absorption_coeff_epidermis(lambda, Cm)
%Calculates absorption coefficient in the epidermis
%with the empirical formula described in the following paper:
% Effects of wavelength-dependent fluence attenuation on the noninvasive
% photoacoustic imaging of hemoglobin oxygen saturation in subcutaneous 
% vasculature in vivo, Inverse Problems 23 (2007) S113-S122.

% INPUT: - lambda, wavelength in [nm]
%        - Cm, volume fraction of melanosome, varies from 1.3% to 43%
%              depending on the skin color, in [%]
% OUTPUT: absorpCoeffEpi, absorption coefficient in the epidermis in [cm-1]

Cm = Cm / 100;
absorpCoeffEpi = (Cm*6.6)*((10^11)*(1/(lambda^3.33))) + ...
                 ((1-Cm)*7.84)*((10^7)*(1/(lambda^3.255)));

end