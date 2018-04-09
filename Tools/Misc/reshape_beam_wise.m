function [outputSignal,uniqueBeamCodes,outputTime] = reshape_beam_wise(inputSignal,beamCode,inputTime)
%reshape_beam_wise.m Reshapes the inputSignal along beam-wise (careful)
% Input:
%   beamCode   : nPulsesPerRecord x nRecords
%   inputSignal: nSamples x nPulsesTotal
% Output:
%   outputSignal:nSamples-1 x nPulsesInBeam x nUniqueBeams (
% ??? need to standardize this the matrix dimensions!

% Identifying the unique beam codes
uniqueBeamCodes = num2cell(unique(beamCode(:)));
% Identifying the beam indices corresponding to each index of inputSignal
beamIndx = cell2mat((cellfun(@(x) beam_filter(beamCode(:),x),uniqueBeamCodes,'UniformOutput',false ))'); 
nSamples=size(inputSignal,1)-1; %removing the last sample to match with SRI output

nPulsesInBeam=size(beamIndx,1);
nBeams=size(beamIndx,2);

outputSignal=reshape(inputSignal(1:nSamples,beamIndx),nSamples,nPulsesInBeam,nBeams);
uniqueBeamCodes=cell2mat(uniqueBeamCodes);
outputTime=reshape(inputTime(beamIndx),1,nPulsesInBeam,nBeams);
% beamCodes=reshape(beamCode(beamIndx),1,nPulsesInBeam,nBeams);

end
function beamIndx=beam_filter(beamCodeArray,beamCode)
    %Finds the index of array corresponding to given beam code
    beamIndx = find(beamCodeArray==beamCode);
end
