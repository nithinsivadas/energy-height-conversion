function [outputSignal,outputTime] = integrate_pulses(inputSignal,inputTime,nPulsesIntegrated)
%integrate_pulses Integrates the pulses of a decoded, and reshaped raw
%voltages for the number of pulses prescribed.
%   Detailed explanation goes here
nPulsesInBeam = size(inputSignal,2);
nRecords = nPulsesInBeam./nPulsesIntegrated;
nBeams = size(inputSignal,3);
nSamples = size(inputSignal,1);

outputSignal = squeeze(sum(reshape(permute(inputSignal,[2 3 1]),...
    nPulsesIntegrated,nRecords,nBeams,nSamples)));
outputTime.min = squeeze(min(reshape(permute(inputTime,[2 3 1]),...
    nPulsesIntegrated,nRecords,nBeams,1)));
outputTime.max = squeeze(max(reshape(permute(inputTime,[2 3 1]),...
    nPulsesIntegrated,nRecords,nBeams,1)));
end