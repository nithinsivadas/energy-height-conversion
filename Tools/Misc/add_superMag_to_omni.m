function [info] = add_superMag_to_omni(superMagFileStr, omniH5Str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    omniH5Str = 'G:\My Drive\Research\Projects\Data\omni.h5';
end

if nargin < 1 || isempty(superMagFileStr)
    superMagFileStr = 'G:\My Drive\Research\Projects\Data\20201217-07-37-supermag.txt';
end
        
% Parse superMag text file
[data,info] = parse_superMag_data(superMagFileStr);
% Convert data gaps to nan values
data.SML(data.SML == 999999) = nan; %999999 is the data gap identifier
data.SMU(data.SMU == 999999) = nan;

% Display commentary from the superMag text file
disp(info);

% Access omni time array
omni.time = h5read(omniH5Str,'/Time');

% Interpolate and identify the SML and SMU values for each omni time
% instant
F.SML = griddedInterpolant(data.time,data.SML,'nearest','none');
F.SMU = griddedInterpolant(data.time,data.SMU,'nearest','none');

SML = F.SML(omni.time);
SMU = F.SMU(omni.time);

% Adding the SML and SMU variable into the omni.h5 data file
try
    
    write_h5_dataset(omniH5Str,'/Indices/SML',SML',1);
    write_h5_dataset(omniH5Str,'/Indices/SMU',SMU',1);
    write_h5_dataset_attribute(omniH5Str,'/Indices/SML','SML - SuperMag Data Revision 5','[nTimex1]','[nT]');
    write_h5_dataset_attribute(omniH5Str,'/Indices/SMU','SMU - SuperMag Data Revision 5','[nTimex1]','[nT]');
    disp(['Successfully wrote SML and SMU data into' omniH5Str]);

catch ME
    getReport(ME);
    throw(ME);
end

end

