%% Generate Calibration HDF5 for DASC

folderStr = 'G:\My Drive\Research\Projects\Paper 3\Data\DASC_Calibration_Files';

sensorloc=[65.1260,-147.4789,689 ];
projectionAlt = 110; %km
minElevationAngle = 15; %deg
%% Before 2011
AZ_2011 = fitsread(strcat(folderStr,filesep,'PKR_20111006_AZ_10deg.FITS'));
EL_2011 = fitsread(strcat(folderStr,filesep,'PKR_20111006_EL_10deg.FITS'));
[dataFlag_2011, lat_2011, lon_2011, alt_2011] = DASC_aer_to_geodetic_v2019(AZ_2011,EL_2011, 15, 110, sensorloc);
AZ_2011 = rotate_array(AZ_2011,-90);
%%
AZ_2011_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_20110112_AZ_10deg.FITS'));
EL_2011_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_20110112_EL_10deg.FITS'));
[dataFlag_2011_2015, lat_2011_2015, lon_2011_2015, alt_2011_2015] = DASC_aer_to_geodetic_v2019(AZ_2011_2015,EL_2011_2015, 15, 110, sensorloc);
AZ_2011_2015 = rotate_array(AZ_2011_2015,-90);
%%
AZ_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_Az.FIT'));
EL_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_El.FIT'));
[dataFlag_2015, lat_2015, lon_2015, alt_2015] = DASC_aer_to_geodetic_v2019(AZ_2015,EL_2015, 15, 110, sensorloc);
%%
time = [1,datenum('01 Jan 2006'),datenum('06 Oct 2011');...
    2,datenum('07 Oct 2011'),datenum('12 Feb 2015');
    3,datenum('12 Feb 2015'),datenum('02 Sep 2019')];

%%
h5File = fullfile(folderStr,'dasc_calibration.h5');
write_h5_dataset(h5File,'/timeRange',time,1,1);
write_h5_dataset(h5File,'/1/AZ',AZ_2011,0);
write_h5_dataset(h5File,'/1/EL',EL_2011,0);
write_h5_dataset(h5File,'/1/lat',lat_2011,0);
write_h5_dataset(h5File,'/1/lon',lon_2011,0);
write_h5_dataset(h5File,'/1/minElFlag',dataFlag_2011,0);

write_h5_dataset(h5File,'/2/AZ',AZ_2011_2015,0);
write_h5_dataset(h5File,'/2/EL',EL_2011_2015,0);
write_h5_dataset(h5File,'/2/lat',lat_2011_2015,0);
write_h5_dataset(h5File,'/2/lon',lon_2011_2015,0);
write_h5_dataset(h5File,'/2/minElFlag',dataFlag_2011_2015,0);

write_h5_dataset(h5File,'/3/AZ',AZ_2015,0);
write_h5_dataset(h5File,'/3/EL',EL_2015,0);
write_h5_dataset(h5File,'/3/lat',lat_2015,0);
write_h5_dataset(h5File,'/3/lon',lon_2015,0);
write_h5_dataset(h5File,'/3/minElFlag',dataFlag_2015,0);

