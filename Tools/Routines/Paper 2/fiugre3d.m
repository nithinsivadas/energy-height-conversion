%% Paper 2: Figure 3d-e: Plot NOAA-17 flux
clear all;
%% Import data
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
% Load POES
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);

%% Calculating L-shells
omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,'26-Mar-2008 08:00','26-Mar-2008 13:00');


%%
smoothFactor = 10;
kext = find_irbem_magFieldModelNo('TS96');
maginput = filter_irbem_maginput(kext,maginput);

% POES trajectory
poesTimeMinStr = '26-Mar-2008 11:27';
poesTimeMaxStr = '26-Mar-2008 11:33';
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
    thisMaginput = interp1(magTime,maginput,poestime(i));
    poesLm(i,1) = onera_desp_lib_make_lstar(kext,[0,0,0,0,0],2,poestime(i),...
        poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),thisMaginput);  
    poesNFoot(i,:)=geopack_find_foot_point(kext,50,2,poestime(i),...
        poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),85,+1,thisMaginput);
end


%% 
cLim = [6 12];
timeMinStr = '26-Mar-2008 11:27';
timeMaxStr = '26-Mar-2008 11:33';

h = figure;
p = create_panels(h,'totalPanelNo',2,'demargin',15);

% L-Shell variations
p(1,1).select();
p1=plot(poestime,poes.mep0E1(poesTimeIndx)-poes.mep0E3(poesTimeIndx),'r','LineWidth',1);
hold on;
p2=plot(poestime,poes.mep90E1(poesTimeIndx)-poes.mep90E3(poesTimeIndx),'k','LineWidth',1);
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
set(gca,'XTickLabel','','YScale','log');
[TTick,TTickLim] = label_time_axis(poestime,false,1/60,timeMinStr,timeMaxStr);
legend([p1 p2], {'Loss-cone','Trapped'});
ylabel({'NOAA-17','30-300 keV e^-'});

p(1,2).select();
plot(poestime,poes.mep0P6(poesTimeIndx),'r','LineWidth',1);
hold on;
plot(poestime,poes.mep90P6(poesTimeIndx),'k','LineWidth',1);
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
set(gca,'XTickLabel','','YScale','log');
[TTick,TTickLim] = label_time_axis(poestime,true,1/60,timeMinStr,timeMaxStr);
ylabel({'NOAA-17','800-2500 keV p^+'});
add_horizontal_axes(TTick,TTickLim,poestime,poesNFoot(:,2), 'NorthFP Lat', 2);
