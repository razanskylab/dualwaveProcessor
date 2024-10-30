function fluenceRate = calculate_fluence_rate(ppe, beamDiameter)
%Calculate fluence rate in unit of [mJ / cm2]

% convert per-pulse energy into [mJ]
ppe = ppe * 1e3;

% convert beam diameter into [cm]
beamDiameter = beamDiameter * 1e2;

beamArea = pi * (beamDiameter/2)^2; % [cm2]
fluenceRate = ppe / beamArea; % [mJ / cm2]

end
