%% Routine Keogram
% One goal:
% 1. Track the temporal evolution of the arcs as they move equatorward

clear all;
opengl('save', 'software');
if ispc
    h5FileStrMCGR = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
    h5FileStrDASC = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3_smaller_time_range.h5';
    omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
else
    h5FileStrMCGR = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-full_v2.h5';
    omniH5FileStr = '/home/nithin/Documents/git-repos/LargeFiles/omni/omni.h5';
end

%%
timeMinStr = '26 Mar 2008 10:30';
timeMaxStr = '26 Mar 2008 11:44';

siteStr = 'gako';
meridianAngle = -21; %Az from true north

dascData = get_2D_plot_inputs_time_independent(h5FileStrMCGR,...
    'plotModeStr','OpticalImage','site','pokerFlat');
gboData = get_2D_plot_inputs_time_independent(h5FileStrMCGR,...
    'plotModeStr','OpticalImage','site',siteStr);

lat = reshape(dascData.latitude,1,[])';
latPixelNum = ceil(sqrt(length(lat)));
dasc_parallels = linspace(min(lat),max(lat),latPixelNum);
dascData.sensorloc(3) = dascData.sensorloc(3)/1000; %Convert m to km
dascData.meridian = get_magnetic_meridian(dascData.sensorloc,...
    datenum('26 Mar 2008 11:00'),dasc_parallels,110,meridianAngle);


lat = reshape(gboData.latitude,1,[])';
latPixelNum = ceil(sqrt(length(lat)));
gbo_parallels = linspace(min(lat),max(lat),latPixelNum);
% gboData.sensorloc(2) = convert_longitude(gboData.sensorloc(2),'360to180');
gboData.sensorloc(3) = gboData.sensorloc(3)/1000; %Convert m to km
% gboData.longitude = convert_longitude(gboData.longitude,'360to180');
gboData.meridian = get_magnetic_meridian(gboData.sensorloc,...
    datenum('26 Mar 2008 11:00'),gbo_parallels,110,meridianAngle);

figure; 
plot(dascData.meridian,dasc_parallels);
hold on;
plot(gboData.meridian,gbo_parallels);
legend('PokerFlat','GAKO');
%% Crop Image
minIndx = find_time(dascData.time,timeMinStr);
maxIndx = find_time(dascData.time,timeMaxStr);
k = 1;
for iTime = minIndx:1:maxIndx
    dascImageCrop(k,:,:) = read_h5_variable_at_time_v2(h5FileStrMCGR,'/DASC/ASI',iTime);
    k = k + 1;
end
dascTime = dascData.time(minIndx:maxIndx);

%%
minIndx = find_time(gboData.time,timeMinStr);
maxIndx = find_time(gboData.time,timeMaxStr);
k = 1;
for iTime = minIndx:1:maxIndx
    gboImageCrop(k,:,:) = read_h5_variable_at_time_v2(h5FileStrMCGR,['/',upper(siteStr),'/ASI'],iTime)';
    k = k + 1;
end
gboTime = gboData.time(minIndx:maxIndx);
%%
[dascKeo,dascLat,dascMeridian] = create_keogram(dascImageCrop,dascData.latitude,dascData.longitude,...
    'time',datenum('26 Mar 2008 11:00'),'meridian',dascData.meridian);

%%
[gboKeo1,gboLat,gboMeridian] = create_keogram(gboImageCrop,gboData.latitude,gboData.longitude,...
    'time',datenum('26 Mar 2008 11:00'),'meridian',gboData.meridian);
%%
match = 0.05;
gboKeo = match*(gboKeo1-gboData.background)+dascData.background-10;
dascLatIndx = find_altitude(dascLat,63):find_altitude(dascLat,68);
gboLatIndx = find_altitude(gboLat,58):find_altitude(gboLat,63);
finalLat = [gboLat(gboLatIndx),dascLat(dascLatIndx)];
finalLon = [gboMeridian(gboLatIndx),dascMeridian(dascLatIndx)];
[~,index] = sortrows(finalLat');
finalParallels = linspace(min(finalLat),max(finalLat),2*1024);
for iTime = 1:1:length(gboTime)
    idascTime = find_time(dascTime,datestr(gboTime(iTime)));
    tempKeo = [gboKeo(gboLatIndx,iTime)',dascKeo(dascLatIndx,idascTime)'];
    F = griddedInterpolant(finalLat(index),tempKeo(index),'nearest');
    finalKeo(:,iTime) = F(finalParallels);
end
finalTime = gboTime;
%% Plot Keogram
totalPanelNo=3;
% clf;
clim = [340, 400];
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',50);
dt = 15*1/60;
q=p(1);
% 4a. Keogram, White light image
q(2).select();
% colormap(gca,get_colormap('k',[0,1,0]));
colormap(viridis);
ax=plot_2D_time_series(gboTime,gboLat,gboKeo,dt,0,timeMinStr,timeMaxStr);
colorbar_thin;
caxis(clim);
ylabel('GAKO');
ylim([58,66]);

q(1).select();
% colormap(gca,get_colormap('k',[0,1,0]));
ax=plot_2D_time_series(dascTime,dascLat,dascKeo,dt,0,timeMinStr,timeMaxStr);
colorbar_thin;
ylabel('DASC');
caxis([340 400]);
ylim([61,68]);

q(3).select();
% colormap(gca,get_colormap('k',[0,1,0]));
ax=plot_2D_time_series(finalTime,finalParallels,finalKeo,dt,0,timeMinStr,timeMaxStr);
colorbar_thin;
ylabel('Combined');
caxis(clim);
ylim([58,68]);
label_time_axis(finalTime, true,dt,timeMinStr,timeMaxStr);


%% Plot the meridians
dascFov = create_dasc_fov(22.5, 110, dascData);
gboFov = create_dasc_fov(22.5, 110, gboData);
timeStr = '26 Mar 2008 11:01:30';
time = datenum(timeStr);
nTime = length(time);
latLim = [58,70];
lonLim = [-170,-120];
deltaLat = 2;
deltaLon = 10;
h=figure('visible','on');
[ax1]=combine_2D_plots_v3(h5FileStrMCGR,h,...
    'maps',{'OpticalImage','OpticalImage'},...
    'sites',{'pokerFlat',siteStr},...
    'thisTime',time,...
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
hold on;
plotm(dascFov.lat,dascFov.lon,'Color',[0.5,0.5,0.5]);
plotm(gboFov.lat,gboFov.lon,'Color',[0.5,0.5,0.5]);
hold on;
plotm(gbo_parallels,gboData.meridian,'r');
plotm(dasc_parallels,dascData.meridian,'w');
%%
function [fov] = create_dasc_fov(elCutOff, projectionAltitude,dascData)

    az = 0:0.1:360;
    el = elCutOff*ones(size(az));
    alt = projectionAltitude;
    C = define_universal_constants;
    slantRange = (-C.RE.*sind(el) + sqrt((C.RE^2).*(sind(el)).^2 + alt.*1000.*(alt.*1000+2.*C.RE)))./1000; 
    
%     slantRange = alt./sind(el);
    [fov.lat,fov.lon,fov.alt] = aer2geodetic(az,el,slantRange,dascData.sensorloc(1),dascData.sensorloc(2),dascData.sensorloc(3)./1000,wgs84Ellipsoid('km'));
    
end