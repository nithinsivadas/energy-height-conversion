
%% THEMIS ASI plotting script
clear all;

%% Initializing
thmasiCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\ThemisASI';

thmasiDataFile = 'thg_l1_asf_gako_2008032611_v01.cdf';
thmasiCalFile = 'thg_l2_asc_gako_19700101_v01.cdf';
locFile = ['C:\Users\nithin\Documents\GitHub\energy-height-conversion\',...
    'Tools\External Tools\thmasi\THEMIS_ASI_Station_List_Nov_2011.xls'];
% [data,dataInfo] = spdfcdfread([thmasiCDFPath,filesep,thmasiDataFile]);
% [cal,calInfo] = spdfcdfread([thmasiCDFPath,filesep,thmasiCalFile]);
fileinfo = parse_thg_filename(thmasiDataFile);
h5OutputFile = 'temp.h5';

%% 
thgdata = parse_thg_cdfData([thmasiCDFPath,filesep,thmasiDataFile],[thmasiCDFPath,filesep,thmasiCalFile]);

site=parse_thg_location_xls(locFile);

siteID=find(strcmpi(site.code,fileinfo.site));
%%
write_thg_to_hdf5(h5OutputFile,thgdata.ASI,'lat',thgdata.glat,...
    'lon',thgdata.glon,'az',thgdata.az,'el',thgdata.el,'alt',thgdata.alt,...
    'altIndx',2,'mlat',thgdata.mlat,'mlon',thgdata.mlon,...
    'time',thgdata.time,'sensorloc',[site.glat(siteID),site.glon(siteID),0],...
    'siteCode',fileinfo.site);
