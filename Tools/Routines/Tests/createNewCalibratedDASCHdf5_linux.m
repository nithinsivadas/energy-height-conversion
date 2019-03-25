%% Testing writing new DASC Calibrated data
clear all;
 localStoreDir = '/home/nithin/Documents/git-repos/Largefiles/';
%  outputH5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\2008_03_26_PKR_DASC.h5';
 outputH5FileStr = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-full_v2.h5';
 projectionAltitudeDASC = 110;
 minElevation = 10;
 minTimeStr = '26 Mar 2008';
 maxTimeStr = [];
% minTimeStr = '26 Mar 2008 10:30:00';
% maxTimeStr = '26 Mar 2008 12:00:00';
calFileLon = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glon.dat'];
calFileLat = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glat.dat'];
setDownloadDASCFlag = false;

%create_DASC_hdf5_low_memory_v2(localStoreDir, outputH5FileStr, projectionAltitudeDASC, minElevation, minTimeStr, maxTimeStr, calFileLat, calFileLon, setDownloadDASCFlag);
%%
add_ASI_background_to_hdf5('dasc',outputH5FileStr); % Adding the backgroud

