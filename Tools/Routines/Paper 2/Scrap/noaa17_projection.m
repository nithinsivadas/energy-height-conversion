%% For Paper 2, Figure 1b; Reviewer 2 - Checking mapping of NOAA-17 using measurements of 0.05 - 20 keV
% Generate figure with satellite tracks, and Themis All Sky Cameras
%% Load Data
% load('G:\My Drive\Research\Projects\Collaborations\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
omnih5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');
magFieldNo = find_irbem_magFieldModelNo('TS96');
maginput = filter_irbem_maginput(magFieldNo,maginput);

%%
% POES trajectory
poesTimeMinStr = '26-Mar-2008 11:27';
poesTimeMaxStr = '26-Mar-2008 11:33';
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
    thisMaginput = interp1(timeMaginput,maginput,poestime(i));
    poesNFoot(i,:)=geopack_find_foot_point(magFieldNo,50,2,poestime(i),...
        poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),85,+1,thisMaginput);
%     poesNFoot(i,:)=onera_desp_lib_find_foot_point(magFieldNo,[0,0,0,0,0],...
%         2,poestime(i),poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),110,+1,thisMaginput);
end
%%
tIndxOffset = -10;
poesColor = poes.mep0E1(poesTimeIndx+tIndxOffset)-poes.mep0E3(poesTimeIndx+tIndxOffset);
poesColorTrapped = poes.mep90E1(poesTimeIndx+tIndxOffset)-poes.mep90E3(poesTimeIndx+tIndxOffset);
poesAnisotropy = poesColor./poesColorTrapped;

poesColor1 = poes.ted01(poesTimeIndx+tIndxOffset);
poesColor2 = poes.ted03(poesTimeIndx+tIndxOffset);
poesColor3 = poes.tedfx5(poesTimeIndx);
poesColor4 = poesColor1.*(1.756*10^-6)+poesColor2.*(5.328*10^-5);
% Find northern footpoint

figure; 
semilogy(datetime(datevec(poes.time(poesTimeIndx))),poesColor4); 
ylim([0.01,10]);
set(gca,'YTick',[0.01,0.1,1,10]);
ylabel({'TED','0.05-20 keV','mW/m^2'});


%% thd
% timeThd = padData.thd.time;
% timeMinStr = '26 Mar 2008 08:00';
% timeMaxStr = '26 Mar 2008 12:00';

% latThd = conv(padData.thd.latFoot,ones(5,1)/5,'same');
% lonThd = conv(padData.thd.lonFoot,ones(5,1)/5,'same');

% latThd = padData.thd.latFoot;
% lonThd = padData.thd.lonFoot;

% latThdL = conv(padData.thd.latFoot,ones(8,1)/8,'same');
% lonThdL = conv(padData.thd.lonFoot,ones(8,1)/8,'same');


% timeMinIndx = find_time(timeThd,timeMinStr);
% timeMaxIndx = find_time(timeThd,timeMaxStr);
% tIndxThd = timeMinIndx:timeMaxIndx;

%% the
% timeThe = padData.the.time;
% timeMinStr = '26 Mar 2008 08:00';
% timeMaxStr = '26 Mar 2008 12:00';

% latThd = conv(padData.the.latFoot,ones(5,1)/5,'same');
% lonThd = conv(padData.the.lonFoot,ones(5,1)/5,'same');

% latThe = padData.the.latFoot;
% lonThe = padData.the.lonFoot;

% latTheL = conv(padData.the.latFoot,ones(8,1)/8,'same');
% lonTheL = conv(padData.the.lonFoot,ones(8,1)/8,'same');


% timeMinIndx = find_time(timeThe,timeMinStr);
% timeMaxIndx = find_time(timeThe,timeMaxStr);
% tIndxThe = timeMinIndx:timeMaxIndx;

%% Initalize
[fov] = create_dasc_fov(fileStr, 22.5, 110);
timeMinStr = '26 Mar 2008 11:29:30';
timeMaxStr = '26 Mar 2008 11:29:30';
time = datenum(timeMinStr):5/(24*60*60):datenum(timeMaxStr);
nTime = length(time);
latLim = [56,72];
lonLim = [-170,-130];
deltaLat = 2;
deltaLon = 10;
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\Figure1b\';
for i = 1:1:nTime
h=figure('visible','on');
[ax1]=combine_2D_plots_v3(fileStr,h,...
    'maps',{'OpticalImage','OpticalImage'},...
    'sites',{'pokerFlat','mcgr'},...
    'thisTime',time(i),...
    'latLim',latLim,...
    'lonLim',lonLim,...
    'elCutOff',5,...
    'deltaLat',deltaLat,...
    'deltaLon',deltaLon,...
    'opticalLim',[0 1],... %[250 450]
    'peakIonizationAltitude',85,...
    'transparency',0.9,...
    'setStoreImage',false);
cm=get_colormap('k',[0.9,1,0.9]);
colormap(viridis);
caxis([0.5,1]);

hold on;
plotm(fov.lat,fov.lon,'Color',[0.5,0.5,0.5]);
% THD
% scatterm(latThd(tIndxThd),lonThd(tIndxThd),6,'b','filled');
% hold on;
% plotm(latThdL(tIndxThd),lonThdL(tIndxThd),'b','LineWidth',1);
% thisTimeStr = datestr(time(i));
% thisTimeIndxThd = find_time(timeThd,thisTimeStr);
% scatterm(latThd(thisTimeIndxThd),lonThd(thisTimeIndxThd),7,'b');
% hhmm = {'08:00','09:00','10:00','10:30'};
% tArrayThd = datenum(strcat('26 Mar 2008'," ",hhmm));
% plot_time_markers(timeThd,latThd,lonThd,tArrayThd,'Northeast','k',[0 0 0]);

% THE
% scatterm(latThe(tIndxThe),lonThe(tIndxThe),6,'c','filled');
% hold on;
% plotm(latTheL(tIndxThe),lonTheL(tIndxThe),'c','LineWidth',1);
% thisTimeStr = datestr(time(i));
% thisTimeIndxThe = find_time(timeThe,thisTimeStr);
% scatterm(latThe(thisTimeIndxThe),lonThe(thisTimeIndxThe),7,'c');
hhmm = {'08:00','09:00','10:00','10:30'};
tArrayThe = datenum(strcat('26 Mar 2008'," ",hhmm));
% plot_time_markers(timeThe,latThe,lonThe,tArrayThe,'Southwest','k',[0 0 0]);

% POES NOAA17
% Anisotropy
ax2=axes;
axm2 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax2.Color = 'none';
ax2.Visible = 'off';
cmap1 = colormap(ax2,get_colormap('w','b'));



% scatterm(poesNFoot(:,2)'-0.15,poesNFoot(:,3)'+0.05...
%     ,10,poesAnisotropy','filled');
scatterm(poesNFoot(:,2)'-0.15,poesNFoot(:,3)'+0.05,...
    10,log10(poesColor'),'filled');

thisPoesNFoot = interp1(poestime,poesNFoot,time(i));

cb2 = colorbar('Location','southoutside');
cb2.Position(2) = 0.2;
cb2.Position(3) = 0.25;
cb2.Position(4) = 0.01;
% cb2.Label.String = {'NOAA17 electrons 30-300 keV','Anisotropy (\phi_|_|/\phi_\perp) [a.u.]'};
cb2.Label.String = {'NOAA17 precipitating electrons 30-300 keV','log10 [Counts/s]'};
% ax2.CLim = [0.5 1];
ax2.CLim = [2.5 3.5];

% Flux
ax3=axes;
axm3 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax3.Color = 'none';
ax3.Visible = 'off';
cmap3 = colormap(ax3,get_colormap('w','r'));

% scatterm(poesNFoot(:,2)',poesNFoot(:,3)',...
%     10,log10(poesColor4'),'filled');

scatterm(poesNFoot(:,2)',poesNFoot(:,3)',...
    10,(poesColor4'),'filled');

hold on;
plotm(thisPoesNFoot(:,2),thisPoesNFoot(:,3),'or');

cb3 = colorbar('Location','southoutside');
cb3.Position(1) = 0.5;
cb3.Position(2) = 0.2;
cb3.Position(3) = 0.25;
cb3.Position(4) = 0.01;
cb3.Label.String = {'NOAA17 precipitating electrons 1-20 keV','[mW/m^2 sr]'};
% ax3.CLim = [2.5 3.5];
% ax3.CLim = [-2 0];
ax3.CLim = [0.1 0.5];

resize_figure(h,200,350);
ax2 = copy_axes_properties(ax1,ax2);
ax3 = copy_axes_properties(ax1,ax3);



imageName = ['figure1b_',datestr(time(i),'HH_MM_SS')];
% export_fig(strcat(storeImageDir,imageName,'_with_POES.png'),'-r600','-png','-nocrop');
end






function plot_time_markers(time,lat,lon,markTimeArr,orientation,tickColor,textColor)
    if nargin<7
        textColor = 'k';
    end
    if nargin <6
        tickColor = 'k';
    end

    if strcmpi(orientation(1:5),'North')
%         dlat = +0.2;
        dlat = +0.05;
    elseif strcmpi(orientation(1:5),'South')
%         dlat = -0.2;
        dlat = -0.05;
    end
    if strcmpi(orientation(6:9),'East')
%         dlon = +0.5;
        dlon = +0.1;
    elseif strcmpi(orientation(6:9),'West')
%         dlon = -3;
        dlon = -1;
    end

    for i = 1:length(markTimeArr)
        thisTimeIndx = find_time(time,datestr(markTimeArr(i)));
        plotm(lat(thisTimeIndx),lon(thisTimeIndx),'.','Color',tickColor);
        textm(lat(thisTimeIndx)+dlat,lon(thisTimeIndx)+dlon...
            ,datestr(markTimeArr(i),'HH:MM'),...
            'Color',textColor);
                %         create_perp_line(lat,lon,thisTimeIndx,orientation(1:5));
    end

end

function create_perp_line(lat,lon,markTimeArrIndx,orientation,tickLength)
    if nargin<6
        tickLength = 0.1;
    end

    if strcmpi(orientation(1:5),'North')
        dlat = 1;
    elseif strcmpi(orientation(1:5),'South')
        dlat = -1;
    end

    for i = length(markTimeArrIndx)
        thisTimeIndx = markTimeArrIndx(i);
        lat1 = lat(thisTimeIndx);
        lon1 = lon(thisTimeIndx);
        dxlon = mean(diff(lon(thisTimeIndx-1:thisTimeIndx+1)));
        dylat = gradient(lat(thisTimeIndx-1:thisTimeIndx+1),dxlon);
        slope = -1./dylat(2);
        c = lat1-slope*lon(1);
        lat2 = lat1+dlat*tickLength;
        lon2 = (lat2-c)./slope;
        hold on;
        plotm([lat1,lat2],[lon1,lon2],'k');
    end

end

function [fov] = create_dasc_fov(fileStr, elCutOff, projectionAltitude)

    dascData = get_2D_plot_inputs_time_independent(fileStr,'plotModeStr','OpticalImage','site','pokerFlat');
    az = 0:0.1:360;
    el = elCutOff*ones(size(az));
    alt = projectionAltitude;
    C = define_universal_constants;
    slantRange = (-C.RE.*sind(el) + sqrt((C.RE^2).*(sind(el)).^2 + alt.*1000.*(alt.*1000+2.*C.RE)))./1000; 
    
%     slantRange = alt./sind(el);
    [fov.lat,fov.lon,fov.alt] = aer2geodetic(az,el,slantRange,dascData.sensorloc(1),dascData.sensorloc(2),dascData.sensorloc(3)./1000,wgs84Ellipsoid('km'));
    
end