%% Routine to test plotting combined 2D plots
%% Initializing
clear all;
homeDir = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\';
inputH5FileStr=[homeDir,'20080326.001_bc_15sec-energyFlux.h5'];
energySlice = 100;
magFieldModelStr = 'TS96';
timeMinStr = '26 Mar 2008 11:05';
timeMaxStr = '26 Mar 2008 11:10';

%% Generating Inputs
pfisrData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','EnergyFluxMap','energySlice',100);

dascData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','OpticalImage');

magData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','MagneticFieldMap','magFieldModelStr',magFieldModelStr);
%% Combine more than one plot for a particular time instant
time = pfisrData.time;
timeMinIndx = find_time(time,timeMinStr);
timeMaxIndx = find_time(time,timeMaxStr);

timeArray = timeMinIndx:1:timeMaxIndx;
% timeArray = 1:1:length(time);

outputImageDir = [homeDir,...
    filesep,'Figure_',...
    num2str(energySlice),'_keV_',datestr(time(timeMinIndx),'yyyymmdd'),filesep];

videoFileNameStr = ['Video_',num2str(energySlice),'_keV_',magFieldModelStr,'_',...
    datestr(time(timeArray(1)),'dd_mmm_yyyy_HH_MM'),'_to_',datestr(time(timeArray(end)),'HH_MM'),'.avi'];

set(0, 'DefaultFigureVisible', 'off');
multiWaitbar('Creating Images',0);
di = 1./length(timeArray);
for itime = timeArray
    try
    h=figure;
    combine_2D_plots(inputH5FileStr,h,...
        'map1','OpticalImage','map3','EnergyFluxMap','map2','MagneticFieldMap',...
    'map1Data',dascData,'map3Data',pfisrData,'map2Data',magData,...
    'energySlice',energySlice,'thisTime',time(itime),'contourLineArray',[1,(10:10:200)],...
    'contourLabelArray',[1,10:10:200],...
    'latLim',[63 67],'lonLim',[-153 -143],'opticalLim',[300 450],...
    'plotContours','Kc','magneticFieldModel',magFieldModelStr,...
    'setStoreImage',true,'imageStoreDir',outputImageDir);
    % 'figureLength',216,'figureBreadth',279);
    catch ME;
    end
    multiWaitbar('Creating Images','Increment',di);
end

set(0, 'DefaultFigureVisible', 'on');

%% Create Video
% tempDir = strsplit(outputImageDir,filesep);
% iImageDir=length(find(~cellfun(@isempty,tempDir)));
% if isunix 
%     iImageDir=iImageDir+1; 
% end
% create_video(strjoin(tempDir(1:iImageDir-1),filesep),tempDir{iImageDir},videoFileNameStr);
% multiWaitbar('Creating Images','Close');