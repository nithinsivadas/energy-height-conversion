%% Download Themis Data

[thmData, matFilePath] = process_themis_data('26-Mar-2008', [initialize_root_path,'LargeFiles',filesep],...
    'tha,thd,the','state');
[maginput,timeMaginput]=generate_maginput('C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5','26 Mar 2008 00:00','27 Mar 2008 00:01');
probes = {'tha','thd','the'};

%% Generate magnetic footprints
timeMinStr = '26-Mar-2008 08:00';
timeMaxStr = '26-Mar-2008 13:00';
magFieldNo = 4;
sysaxes = 1;
nProbes = length(probes);
iP = 1./nProbes;
multiWaitbar('1.Probes...',0);
for iProbe = 1:1:nProbes
    timeIndx = find_time(thmData.(probes{iProbe}).state.time,datenum(timeMinStr)):1:find_time(thmData.(probes{iProbe}).state.time,datenum(timeMaxStr));
    time = thmData.(probes{iProbe}).state.time(timeIndx);
    multiWaitbar('2.Footprint...',0);
    nTime = length(time);
    dt = 1./nTime;
    thmFoot.(probes{iProbe}).magFieldModelStr=find_irbem_magFieldModelStr(magFieldNo);
    thmFoot.(probes{iProbe}).time = time;
    for iTime = 1:1:nTime
        thisMaginput = interp1(timeMaginput',maginput,time(iTime));
        thmFoot.(probes{iProbe}).GDZ(iTime,:) = onera_desp_lib_find_foot_point...
            (magFieldNo,[0,0,0,0,0],sysaxes,time(iTime),...
            thmData.(probes{iProbe}).state.XYZ_GEO(timeIndx(iTime),1),...
            thmData.(probes{iProbe}).state.XYZ_GEO(timeIndx(iTime),2),...
            thmData.(probes{iProbe}).state.XYZ_GEO(timeIndx(iTime),3),...
            110,+1,thisMaginput); 
    multiWaitbar('2.Footprint...','Increment',dt);
    end
    multiWaitbar('1.Probes...','Increment',iP);
end

%% Load NOAA Data
poesCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17';
poesbinName = 'n1720080326.mat';
load([poesCDFPath,filesep,poesbinName]);
noaa17 = poes;
time = noaa17.time;
multiWaitbar('3.NOAA17 Footprint...',0);
nTime = length(time);
dt = 1./nTime;
noaa17Foot.magFieldModelStr=find_irbem_magFieldModelStr(magFieldNo);
noaa17Foot.time = time;
sysaxes = 0;
for iTime = 1:1:nTime
    thisMaginput = interp1(timeMaginput',maginput,time(iTime));
    noaa17Foot.GDZ(iTime,:) = onera_desp_lib_find_foot_point...
        (magFieldNo,[0,0,0,0,0],sysaxes,time(iTime),...
        noaa17.lat(iTime,1),...
        noaa17.lon(iTime,1),...
        850,... %altitude in km
        110,+1,thisMaginput); 
multiWaitbar('3.NOAA17 Footprint...','Increment',dt);
end

%% Plot
inputH5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20080326.001_bc_15sec-energyFlux_v85.h5';
h = figure;
timeStr = '26-Mar-2008 10:00:00';
energySlice = 100;
tic
p = combine_2D_plots_v2(inputH5FileStr,h,...
        'maps',{'OpticalImage','OpticalImage','OpticalImage','EnergyFluxMap','OpticalImage','OpticalImage'},...
    'sites',{'gako','pokerFlat','mcgr','pokerFlat','inuv','whit'},...
    'energySlice',energySlice,...
    'thisTime',datenum(timeStr),...
    'latLim',[58 72],...
    'deltaLat',2,...
    'deltaLon',10,...
    'lonLim',[-160 -125],...
    'opticalLim',[250 350],...
    'peakIonizationAltitude',85,...
    'elCutoff',20,...
    'transparency',0.8,...
    'setStoreImage',false);
toc

hold on;
timeIndx = find_time(thmFoot.thd.time,datenum('26-Mar-2008 08:00')):1:find_time(thmFoot.thd.time,datenum('26-Mar-2008 12:00'));
thisTimeIndx = find_time(thmFoot.thd.time,timeStr);
scatterm(thmFoot.thd.GDZ(timeIndx,2),thmFoot.thd.GDZ(timeIndx,3),3,'r');
scatterm(thmFoot.thd.GDZ(thisTimeIndx,2),thmFoot.thd.GDZ(thisTimeIndx,3),15,'k');
textm(thmFoot.thd.GDZ(thisTimeIndx,2),thmFoot.thd.GDZ(thisTimeIndx,3)-0.3,'thd','HorizontalAlignment','right');

% pfisrData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
%     'plotModeStr','EnergyFluxMap','energySlice',100);
% dascData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
%     'plotModeStr','OpticalImage');
% 
% p = combine_2D_plots(inputH5FileStr,h,...
%         'map1','OpticalImage','map2','EnergyFluxMap',...
%     'map1Data',dascData,'map2Data',pfisrData,...
%     'energySlice',energySlice,...
%     'thisTime',datenum(timeStr),...
%     'latLim',[63 67],...
%     'lonLim',[-153 -143],...
%     'opticalLim',[300 450],...
%     'setStoreImage',false);

% %% Initializing
% thmasiCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\ThemisASI';
% 
% thmasiDataFile = 'thg_l1_asf_gako_2008032611_v01.cdf';
% thmasiCalFile = 'thg_l2_asc_gako_19700101_v01.cdf';
% locFile = ['C:\Users\nithin\Documents\GitHub\energy-height-conversion\',...
%     'Tools\External Tools\thmasi\THEMIS_ASI_Station_List_Nov_2011.xls'];
% % [data,dataInfo] = spdfcdfread([thmasiCDFPath,filesep,thmasiDataFile]);
% % [cal,calInfo] = spdfcdfread([thmasiCDFPath,filesep,thmasiCalFile]);
% fileinfo = parse_thg_filename(thmasiDataFile);
% h5OutputFile = 'temp.h5';
% 
% %% 
% thgdata = parse_thg_cdfData([thmasiCDFPath,filesep,thmasiDataFile],[thmasiCDFPath,filesep,thmasiCalFile]);
% 
% site=parse_thg_location_xls(locFile);
% 
% siteID=find(strcmpi(site.code,fileinfo.site));
% %%
% write_thg_to_hdf5(h5OutputFile,thgdata.ASI,'lat',thgdata.glat,...
%     'lon',thgdata.glon,'az',thgdata.az,'el',thgdata.el,'alt',thgdata.alt,...
%     'altIndx',2,'mlat',thgdata.mlat,'mlon',thgdata.mlon,...
%     'time',thgdata.time,'sensorloc',[site.glat(siteID),site.glon(siteID),0],...
%     'siteCode',fileinfo.site);
