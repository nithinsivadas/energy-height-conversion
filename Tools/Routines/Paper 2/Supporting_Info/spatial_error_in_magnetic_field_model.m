%% NOAA-17 mapping error
%% Load Data
% load('G:\My Drive\Research\Projects\Collaborations\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
omnih5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');

%%
magFieldNo = find_irbem_magFieldModelNo('TS01');
maginput = filter_irbem_maginput(magFieldNo,maginput);

%%
% POES trajectory
timeMinStr = '26 Mar 2008 11:29:30';
timeMaxStr = '26 Mar 2008 11:29:30';
poesTimeMinStr = timeMinStr;
poesTimeMaxStr = timeMaxStr;
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
    thisMaginput = interp1(timeMaginput,maginput,poestime(i));
    poesNFoot(i,:)=geopack_find_foot_point(magFieldNo,50,2,poestime(i),...
        poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),85,+1,thisMaginput);
    kext1 = 1:1:11;
    for k=1:1:length(kext1)
        maginput1 = filter_irbem_maginput(kext1(k),maginput);
        thisMaginput = interp1(timeMaginput,maginput1,poestime(i));
        [~,poesEqPoint1(k,:)] = onera_desp_lib_find_magequator(kext1(k),[0,0,0,0,0],2,poestime(i),...
            poesGSM(i,1),poesGSM(i,2),poesGSM(i,3),thisMaginput);
        poesNFoot1(k,:)=onera_desp_lib_find_foot_point(kext1(k),[0,0,0,0,0],...
            1,poestime(i),poesEqPoint1(i,1),poesEqPoint1(i,2),poesEqPoint1(i,3),85,+1,thisMaginput);
    end
end

marker = ['o','+','*','.','x','s','d','^','v','>','<'];

tIndxOffset = 0;
[fov] = create_dasc_fov(fileStr, 22.5, 110);
time = datenum(timeMinStr):5/(24*60*60):datenum(timeMaxStr);
nTime = length(time);
latLim = [58,68];
lonLim = [-170,-135];
deltaLat = 2;
deltaLon = 10;
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\Figure1b\';
for i = 1:1:nTime
h=figure('visible','on');
[ax1]=combine_2D_plots_v3(fileStr,h,...
    'maps',{'OpticalImage'},...
    'sites',{'mcgr'},...
    'thisTime',time(i),...
    'latLim',latLim,...
    'lonLim',lonLim,...
    'elCutOff',15,...
    'deltaLat',deltaLat,...
    'deltaLon',deltaLon,...
    'opticalLim',[0 1],... %[250 450]
    'peakIonizationAltitude',85,...
    'transparency',0.9,...
    'setStoreImage',false);
cm=get_colormap('k',[0.9,1,0.9]);
colormap(viridis);
caxis([0.5,1]);

% hold on;
% plotm(fov.lat,fov.lon,'Color',[0.5,0.5,0.5]);

% POES NOAA17
% Anisotropy
ax2=axes;
axm2 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax2.Color = 'none';
ax2.Visible = 'off';
for k=1:1:11
    scatterm(poesNFoot1(k,2)',poesNFoot1(k,3)',marker(k),'m');
    hold on;
    scatterm(68-0.5*k,-140.5,marker(k),'k');
    hold on;
    textm(68-0.5*k,-140,find_irbem_magFieldModelStr(k),'Color','k');
end
hold on;
% scatterm(poesNFoot(1,2)',poesNFoot(1,3)','p','r');
% scatterm(68-0.5*12,-140.5,'p','k');
% textm(68-0.5*12,-140,[find_irbem_magFieldModelStr(magFieldNo),' GEOPACK'],'Color','k');
% thisPoesNFoot = interp1(poestime,poesNFoot,time(i)-tIndxOffset*2/(3600*24));
textm(68-0.5*12+0.15,-141.7,'\DeltaLat: ','Color','k');
textm(68-0.5*12,-140,[num2str(max(poesNFoot1(:,2))-min(poesNFoot1(:,2)),'%2.1f'),'°N'],'Color','k');
% cb2 = colorbar('Location','southoutside');
% cb2.Position(2) = 0.2;
% cb2.Position(3) = 0.25;
% cb2.Position(4) = 0.01;
% cb2.Label.String = {'NOAA17 electrons 30-300 keV','Anisotropy (\phi_|_|/\phi_\perp) [a.u.]'};
% ax2.CLim = [0.5 1];

resize_figure(h,200,350);
ax2 = copy_axes_properties(ax1,ax2);



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