%% Testing phase
clear all;

%% Extracting data
poesCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17';
poesCDFName = 'poes_n17_20080326.cdf';
poesnCDFName = 'POES_combinedSpectrum_n17_90_20080326.nc';
poesbinName = 'n1720080326.mat';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';

%% Opening netcdf
load([poesCDFPath,filesep,poesbinName]);
% ncid = netcdf.open([poesCDFPath,filesep,poesnCDFName]);
noaa17 = poes;

%%
[data,info] = cdfread([poesCDFPath,filesep,poesCDFName]);

poesData.time = cellfun(@(x) todatenum(x),data(:,1));
tempgeoLL = cell2mat(data(:,2)');
poesData.lat = tempgeoLL(1,:)';
poesData.lon = tempgeoLL(2,:)';
poesData.Lm = cell2mat(data(:,4));
poesData.MLT = cell2mat(data(:,5));
poesData.totalEnergyFlux = cell2mat(data(:,8));
poesData.mep0P = cell2mat(data(:,12)')'; 
poesData.mep90P = cell2mat(data(:,13)')'; 
poesData.mep0E = cell2mat(data(:,14)')';
poesData.mep90E = cell2mat(data(:,15)')'; 
poesData.eLabel = data{1,18};
poesData.pLabel = data{1,19};

% [maginput,magtime] = generate_maginput(omniH5FileStr, poesData.time(1), poesData.time(end));
% poesData.maginput = interp_nans(interp1(magtime',maginput,poesData.time));

%%
% [maginput1,magtime1] = generate_maginput(omniH5FileStr, noaa17.time(1), noaa17.time(end));
% noaa17.maginput = interp_nans(interp1(magtime',maginput,noaa17.time));

%%
% magFieldModel = 3;
% ntime = find_time(poesData.time,'26-Mar-2008 11:32');
% tic
% [Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar(magFieldModel,[1,0,7,7,0],0,poesData.time(ntime,1),802*ones(size(poesData.lat(ntime,1))),poesData.lat(ntime,1),poesData.lon(ntime,1),poesData.maginput(ntime,:));
% [Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar(magFieldModel,[0,0,0,0,0],0,poesData.time,824*ones(size(poesData.lat)),poesData.lat,poesData.lon,poesData.maginput);
% [Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar(magFieldModel,[0,0,0,0,0],0,poesData.time,850*ones(size(noaa17.lat)),noaa17.lat,noaa17.lon,noaa17.maginput);
% [Xfoot, Bfoot, Bfootmag] = onera_desp_lib_find_foot_point(magFieldModel,[0,0,0,0,0],0,noaa17.time,850*ones(size(noaa17.lat)),noaa17.lat,noaa17.lon,110,+1,noaa17.maginput);
% [Xfoot, Bfoot, Bfootmag] = onera_desp_lib_find_foot_point(magFieldModel,[1,0,7,7,0],0,poesData.time(ntime,1),802*ones(size(poesData.lat(ntime,1))),poesData.lat(ntime,1),poesData.lon(ntime,1),110,+1,poesData.maginput(ntime,:));
% poesData.maginput(ntime,[5,2,6,7])
% toc

%% AAGCMv2
% ver=pyversion;
% multiWaitbar('AACGMV2...',0);
% di = 1./length(poes.time);
% if ver=='3.6'
%     for i=1:1:length(poes.time)
%     multiWaitbar('AACGMV2...','Increment',di);
%     tup=py.aacgmv2.wrapper.get_aacgm_coord(noaa17.lat(i), noaa17.lon(i), 850, noaa17.time(i), 'TRACE');
%     mlat(i) = double(py.array.array('d',py.numpy.nditer(tup{1})));
%     mlon(i) = double(py.array.array('d',py.numpy.nditer(tup{2})));
%     mlt(i) = double(py.array.array('d',py.numpy.nditer(tup{3})));
%     end
% else
%     error('Python not integrated to matlab');
% end

%%
% figure;
% plot(datetime(datestr(poesData.time)),abs(Lstar),'.-'); hold on; 
% plot(datetime(datestr(poesData.time)),poesData.Lm);
% xlim([datetime('26-Mar-2008 11:00'),datetime('26-Mar-2008 11:30')]);
% legend([find_irbem_magFieldModelStr(magFieldModel),' Lstar'],'POES Lm');
% tic
% [Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_make_lstar(3,[0,0,0,0,0],0,poesData.time,824*ones(size(poesData.lat)),poesData.lat,poesData.lon,poesData.maginput);
% toc
%%
% timeMinStr = '26-Mar-2008 11:00';
% timeMaxStr = '26-Mar-2008 12:00';
% timeIndx = find_time(poesData.time,timeMinStr):1:find_time(poesData.time,timeMaxStr);
% figure('Color','w'); axesm('eqaazim','MapLatLimit',[30 90],'MapLonLimit',[-160 -135]); axis off; framem on; gridm on; mlabel on; plabel on; setm(gca, 'MLabelParallel',0); plotm(Xfoot(timeIndx,2),Xfoot(timeIndx,3));
% geoshow(coastlat,coastlon,'DisplayType','polygon');

%%
% magFieldModel =4;
% options = [0,0,0,0,0];
% sysaxes = 0;
% timeID = find_time(poesData.time,'26-Mar-2008 11:40');
% [Lm,Blocal,Bmin,J,POSIT] = onera_desp_lib_trace_field_line(magFieldModel,options,sysaxes,poesData.time(timeID),824,poesData.lat(timeID,1),poesData.lon(timeID,1),poesData.maginput(timeID,:),1); 

%% Initializing
timeMinStr = '26-Mar-2008 11:29';
timeMaxStr = '26-Mar-2008 11:31';
timeIndx = find_time(poesData.time,timeMinStr):1:find_time(poesData.time,timeMaxStr);
timeIndxNOAA = find_time(noaa17.time,timeMinStr):1:find_time(noaa17.time,timeMaxStr);
timeTick = 8*1/3600;



%% Plotting
totalPanelNo = 5;

p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',10,'panelSize',30);

q=p(1);

% Panel 1: Total Energy Flux
q(1).select();
plot(poesData.time(timeIndx),poesData.totalEnergyFlux(timeIndx));
set(gca,'YScale','log');
set(gca,'YLim',[10^-5 10^1]);
ylabel('Total energy flux');
label_time_axis(poesData.time(timeIndx), false, timeTick);

% Panel 2: Electron Flux 30 keV
q(2).select();
% plot(poesData.time(timeIndx),poesData.mep0E(timeIndx,1),'.-r');
% hold on;
% plot(poesData.time(timeIndx),poesData.mep90E(timeIndx,1),'.-k');
% hold on;
plot(noaa17.time(timeIndxNOAA),noaa17.mep0E1(timeIndxNOAA),'-r');
hold on;
plot(noaa17.time(timeIndxNOAA),noaa17.mep90E1(timeIndxNOAA),'-k');
set(gca,'YScale','log');
set(gca,'YLim',[10^0 10^4]);
ylabel('>30 keV e^-');
legend('loss-cone', 'trapped','Location','northwest');
% label_time_axis(poesData.time(timeIndx), false, timeTick);
label_time_axis(noaa17.time(timeIndxNOAA), false, timeTick);

% Panel 3: Electron Flux 100 keV
q(3).select();
% plot(poesData.time(timeIndx),poesData.mep0E(timeIndx,2),'r');
% hold on;
% plot(poesData.time(timeIndx),poesData.mep90E(timeIndx,2),'k');
plot(noaa17.time(timeIndxNOAA),noaa17.mep0E2(timeIndxNOAA),'-r');
hold on;
plot(noaa17.time(timeIndxNOAA),noaa17.mep90E2(timeIndxNOAA),'-k');
set(gca,'YScale','log');
set(gca,'YLim',[10^0 10^4]);
ylabel('>100 keV e^-');
legend('loss-cone', 'trapped','Location','northwest');
% label_time_axis(poesData.time(timeIndx), false, timeTick);
label_time_axis(noaa17.time(timeIndxNOAA), false, timeTick);

% Panel 4: Electron Flux 300 keV
q(4).select();
% plot(poesData.time(timeIndx),poesData.mep0E(timeIndx,3),'r');
% hold on;
% plot(poesData.time(timeIndx),poesData.mep90E(timeIndx,3),'k');
plot(noaa17.time(timeIndxNOAA),noaa17.mep0E3(timeIndxNOAA),'-r');
hold on;
plot(noaa17.time(timeIndxNOAA),noaa17.mep90E3(timeIndxNOAA),'-k');
set(gca,'YScale','log');
set(gca,'YLim',[10^0 10^4]);
ylabel('>300 keV e^-');
legend('loss-cone', 'trapped','Location','northwest');
% label_time_axis(poesData.time(timeIndx), false, timeTick);
label_time_axis(noaa17.time(timeIndxNOAA), false, timeTick);

% Panel 5: Proton Flux 7500 keV
q(5).select();
% plot(poesData.time(timeIndx),poesData.mep0P(timeIndx,6),'r');
% hold on;
% plot(poesData.time(timeIndx),poesData.mep90P(timeIndx,6),'k');
plot(noaa17.time(timeIndxNOAA),noaa17.mep0P6(timeIndxNOAA),'-r');
hold on;
plot(noaa17.time(timeIndxNOAA),noaa17.mep90P6(timeIndxNOAA),'-k');
set(gca,'YScale','log');
set(gca,'YLim',[10^0 10^2]);
ylabel('>6900 keV p^+');
legend('loss-cone', 'trapped','Location','northwest');

% [TTick, TTickLim] = label_time_axis(poesData.time(timeIndx), true, timeTick);
% add_horizontal_axes(TTick,TTickLim,poesData.time(timeIndx),poesData.Lm(timeIndx), 'Lm', 2);
% add_horizontal_axes(TTick,TTickLim,poesData.time(timeIndx),poesData.lat(timeIndx), 'Lat', 3);
% add_horizontal_axes(TTick,TTickLim,poesData.time(timeIndx),Xfoot(timeIndx,2),['NFP_Lat [',find_irbem_magFieldModelStr(magFieldModel),']'], 4);
[TTick, TTickLim] = label_time_axis(noaa17.time(timeIndxNOAA), false, timeTick);
add_horizontal_axes(TTick,TTickLim,noaa17.time(timeIndxNOAA),(mod(noaa17.time(timeIndxNOAA),1)*24*60*60)-(mod(noaa17.time(timeIndxNOAA(1)),1)*24*60*60), 'Time', 1);
add_horizontal_axes(TTick,TTickLim,noaa17.time(timeIndxNOAA),noaa17.mlat(timeIndxNOAA), 'MLat', 2);
add_horizontal_axes(TTick,TTickLim,noaa17.time(timeIndxNOAA),noaa17.mlt(timeIndxNOAA), 'MLT', 3);
% add_horizontal_axes(TTick,TTickLim,noaa17.time(timeIndxNOAA),Xfoot(timeIndxNOAA,2),['NFP_Lat [',find_irbem_magFieldModelStr(magFieldModel),']'], 3);
