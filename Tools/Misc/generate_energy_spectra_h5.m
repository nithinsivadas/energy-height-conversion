function [data,status] = generate_energy_spectra_h5...
    (inputAMISRExpFile,filePath,projectionAlt,energyBin,altLim,outputH5FileStr,...
    timeMinStr,timeMaxStr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 8 || isempty(timeMaxStr)
    timeMaxStr = [];
end
if nargin < 7 || isempty(timeMinStr)
    timeMinStr = [];
end
if nargin < 6 || isempty(outputH5FileStr)
    outputH5FileStr = [filePath,inputAMISRExpFile(1:20),'-energyFlux.h5'];
end
if nargin < 5
    altLim = [60 200];
end

if nargin < 4
    energyBin = logspace(3,6,30);
end

if nargin < 3
    projectionAlt = 70;
end

fileNameStr = [filePath,inputAMISRExpFile];
amisrData = read_amisr(fileNameStr);
amisrData = aer_to_field_aligned_coords(amisrData,projectionAlt);
amisrData = interpolate_to_field_aligned_coords(amisrData,timeMinStr,timeMaxStr);
[dataInv, magcoords, dataInputInv] = get_2D_energy_spectra(amisrData,energyBin',...
   timeMinStr,timeMaxStr,altLim,'magnetic');
data = create_energyFlux_hdf5(dataInv,dataInputInv,outputH5FileStr,'dataInv');
status = ['Successfully stored the file in ',outputH5FileStr];
end

