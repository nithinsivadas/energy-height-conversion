function [kpap,kpapHeader] = get_kpap(yearStr,setForceDownload,localStoreDir,url)
%get_kpap Get KP, AP, and F10.7 Index from geophys parameters. Can be used
%         for msis
% Instead of this you can use the external tool f107_aph.m
%   Detailed explanation goes here
if nargin<4
    url='http://45.79.106.44/geophys_params/';
end
if nargin<3
    localStoreDir = [initialize_root_path,'LargeFiles',filesep,'KP_AP',filesep];
end
if nargin<2
    setForceDownload=false;
end


if ~isfile([localStoreDir,yearStr]) || setForceDownload
    remoteFile=[url,yearStr];
    localFilePath=websave([localStoreDir,yearStr],remoteFile);
else
    localFilePath=[localStoreDir,yearStr];
    warning(['File ',yearStr,' alread exists in ',localStoreDir]);
end

fileID=fopen(localFilePath); 
fileFormat='%2d%2d%2d%4d%2d%2d%2d%2d%2d%2d%2d%2d%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%2d%3.1f%1d%3c%5.1f%1d'; 
kpap = textscan(fileID,fileFormat);
fclose(fileID);
kpapHeader = ...
    ['Year';'Month';'Day';'Bartels Solar Rotation Number';...
    'Number of Days within the Bartels 27-day cycle';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 0000 - 0300 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 0300 - 0600 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 0600 - 0900 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 0900 - 1200 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 1200 - 1500 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 1500 - 1800 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 1800 - 2100 UT.';...
    'Kp or PLANETARY 3-HOUR RANGE INDEX for 2100 - 2400 UT.';...
    'SUM of the eight Kp indices for the day expressed to the nearest third of a unit.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 0000 - 0300 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 0300 - 0600 UT.';...
    'PLANETARY EQUIVALENT AMPLITUDE for 0600 - 0900 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 0900 - 1200 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 1200 - 1500 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 1500 - 1800 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 1800 - 2100 UT.';...
    'ap or PLANETARY EQUIVALENT AMPLITUDE for 2100 - 2400 UT.';...
    'Ap or PLANETARY EQUIVALENT DAILY AMPLITUDE--the arithmetic mean of the days eight ap values.';...
    'Cp or PLANETARY DAILY CHARACTER FIGURE--a qualitative estimate of overall level of magnetic activity for the day determined from the sum of the eight ap amplitudes.  Cp ranges, in steps of one-tenth, from 0 (quiet) to 2.5 (highly disturbed).';...
    'C9--a conversion of the 0-to-2.5 range of the Cp index to one digit between 0 and 9.';...
    'INTERNATIONAL SUNSPOT NUMBER.  Records contain the Zurich number through December 31, 1980, and the International Brussels number thereafter.';...
    'OTTAWA 10.7-CM SOLAR RADIO FLUX ADJUSTED TO 1 AU--measured at 1700 UT daily and expressed in units of 10 to the -22 Watts/meter sq/hertz.  Observations began on February 14, 1947. From that date through December 31, 1973, the fluxes given here do not reflect the revisions Ottawa made in 1966. NOTE: If a solar radio burst is in progress during the observation the pre-noon or afternoon value is used (as indicated by a flux qualifier value of 1 in column 71.';...
    'FLUX QUALIFIER.  "0" indicates flux required no adjustment; "1" indicates flux required adjustment for burst in progress at time of measurement; "2" indicates a flux approximated by either interpolation or extrapolation; and "3" indicates no observation.'...
];


end

