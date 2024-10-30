% File: counts_to_voltage.m
% Weiye Li - Razansky Lab - University and ETH Zurich
% E-mail: weiye.li@uzh.ch
% Date: August 18, 2022

function [signal] = counts_to_voltage(signal, measRange, resolution)
  % converts the signal valued in discrete counts to voltages (0 counts are also 0 voltage)

  % INPUT: signal before unit conversion, measurement range of DAQ in [mV], resolution of DAQ (e.g. 16 for 16-bit DAQ)
  % OUTPUT: signal after unit conversion, now in [mV]

  signal = single(signal);
  lsbSize = measRange / (2^(resolution-1));
  signal = signal .* lsbSize;

end
