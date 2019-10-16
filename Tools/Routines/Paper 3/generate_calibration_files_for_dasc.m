%% Generate Calibration HDF5 for DASC
clear all;

folderStr = 'G:\My Drive\Research\Projects\Paper 3\Data\DASC_Calibration_Files';

sensorloc=[65.1260,-147.4789,689 ];
projectionAlt = 110; %km
minElevationAngle = 15; %deg

input = {1,datenum('01 Jan 2006'),datenum('06 Mar 2008 23:59'),'PKR_20111006_AZ_10deg.FITS', 'PKR_20111006_EL_10deg.FITS';...
         2,datenum('07 Mar 2008'),datenum('27 Mar 2008 23:59'),'PKR_DASC_20111006_AZ_10deg_Nithin.FITS', 'PKR_DASC_20111006_EL_10deg_Nithin.FITS';...
%          3,datenum('28 Mar 2008'),datenum('06 Oct 2011 23:59'),'PKR_20111006_AZ_10deg.FITS', 'PKR_20111006_EL_10deg.FITS';...
         3,datenum('28 Mar 2008'),datenum('12 Feb 2013'), 'PKR_DASC_20110112_AZ_10deg_Nithin.FITS', 'PKR_DASC_20110112_EL_10deg_Nithin.FITS';...
            4,datenum('12 Feb 2013'),datenum('02 Sep 2019'), 'PKR_DASC_0558_20150213_Az.FIT', 'PKR_DASC_0558_20150213_El.FIT'};
%%

for i = 1:1:size(input,1)
    AZ{i} = fitsread(strcat(folderStr,filesep,input{i,4}));
    EL{i} = fitsread(strcat(folderStr,filesep,input{i,5}));
    [dataFlag{i}, lat{i}, lon{i}, alt{i}] = DASC_aer_to_geodetic_v2019(AZ{i},EL{i}, minElevationAngle, projectionAlt, sensorloc);
    time(i,:) = cell2mat(input(i,1:3));
end
        
        
%% Before 2011
% AZ_2011 = fitsread(strcat(folderStr,filesep,'PKR_20111006_AZ_10deg.FITS'));
% EL_2011 = fitsread(strcat(folderStr,filesep,'PKR_20111006_EL_10deg.FITS'));
% [dataFlag_2011, lat_2011, lon_2011, alt_2011] = DASC_aer_to_geodetic_v2019(AZ_2011,EL_2011, 15, 110, sensorloc);
% AZ_2011(AZ_2011==0)=nan;
% AZ_2011 = fliplr(rotate_array(AZ_2011,-63));
%%
% AZ_2011_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_20110112_AZ_10deg_Nithin.FITS'));
% EL_2011_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_20110112_EL_10deg_Nithin.FITS'));
% [dataFlag_2011_2015, lat_2011_2015, lon_2011_2015, alt_2011_2015] = DASC_aer_to_geodetic_v2019(AZ_2011_2015,EL_2011_2015, 15, 110, sensorloc);
% AZ_2011_2015(AZ_2011_2015==0)=nan;
% AZ_2011_2015 = fliplr(rotate_array(AZ_2011_2015,+27));
%%
% AZ_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_Az.FIT')); % Seems correct
% EL_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_El.FIT'));
% [dataFlag_2015, lat_2015, lon_2015, alt_2015] = DASC_aer_to_geodetic_v2019(AZ_2015,EL_2015, 15, 110, sensorloc);
% AZ_2015(AZ_2015==0)=nan;
% AZ_2015 = fliplr(rotate_array(AZ_2015,-63));
%%


%%
h5File = fullfile(folderStr,'dasc_calibration.h5');
write_h5_dataset(h5File,'/timeRange',time,1,1);

for i = 1:1:size(input,1)
    indxStr = num2str(i);
    write_h5_dataset(h5File,['/',indxStr,'/AZ'],AZ{i},0);
    write_h5_dataset(h5File,['/',indxStr,'/EL'],EL{i},0);
    write_h5_dataset(h5File,['/',indxStr,'/lat'],lat{i},0);
    write_h5_dataset(h5File,['/',indxStr,'/lon'],lon{i},0);
    write_h5_dataset(h5File,['/',indxStr,'/minElFlag'],dataFlag{i},0);
end

