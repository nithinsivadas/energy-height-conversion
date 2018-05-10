function [decodedSignal] = decoding(rawVoltages,code)
% decoding.m
%   Decoding a phase coded signal by matchfiltering it with the known phase-coded signal
%   Input:
%   rawVoltages - [nSamples x nRanges x nPulsesPerRecord x nRecords]
%   code        - phase code
if nargin<2
   % Default 13-baud barker code, oversampled by a factor of 2
    code=...
    [+1 +1 +1 +1 +1 +1 +1 +1 +1 +1 -1 -1 -1 -1 +1 +1 +1 +1 -1 -1 +1 +1 -1 -1 +1 +1]; %Oversampled by a factor of 2
end
nRecords = size(rawVoltages,4);
nPulsesPerRecord = size(rawVoltages,3);
nRanges = size(rawVoltages,2);
% Converting the Raw Data matrix into a cell array in order matchfilter without for loops
signal = mat2cell(reshape(squeeze(rawVoltages(1,:,:,:)+1i*rawVoltages(2,:,:,:)),nRanges,nPulsesPerRecord*nRecords),nRanges,ones(1,nPulsesPerRecord*nRecords));
% Barker-code Decoding: Matchfiltering with the oversampled barker code, and converting the cell arrays back to matrices
decodedSignal = cell2mat(cellfun(@(x) conv(x,fliplr(conj(code)),'valid'),signal,'UniformOutput',false ));
end
