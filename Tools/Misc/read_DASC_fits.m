function [imageData,time] = read_DASC_fits(fileNameStr,imageSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Read FITS data from file
dataOldRes = fitsread(fileNameStr);
% initialize_geodata;
%% Change resolution of data location and data to size of az & el
imageData = modify_matrix_size(dataOldRes, imageSize, imageSize);

%% Generating Time Stamp
splitFileNameStr = strsplit(fileNameStr,filesep);
time = fitsfiletimestamp(splitFileNameStr(end));
end

