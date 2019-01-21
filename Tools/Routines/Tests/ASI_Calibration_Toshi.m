%% Toshi's Calibration of DASC, Poker Flat
clear all;
%%
latCalFile = 'G:\My Drive\Research\Projects\Data\DASC_PKR_Calibration\200803261040uafdasc_pkr_glat.dat';
lonCalFile = 'G:\My Drive\Research\Projects\Data\DASC_PKR_Calibration\200803261040uafdasc_pkr_glon.dat';

% latCal = dlmread(latCalFile);
% lonCal = dlmread(lonCalFile);

fileLatID = fopen(latCalFile);
latCal = cell2mat(textscan(fileLatID,repmat('%8.4f',1,1024)));
latCal=interp_nans(latCal)';
fclose(fileLatID);

fileLonID = fopen(lonCalFile);
lonCal = cell2mat(textscan(fileLonID,repmat('%8.4f',1,1024)));
lonCal=interp_nans(lonCal)';
fclose(fileLonID);

%%
stars=get_star_catalogue;
% fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\testDASC.h5';

%% Step 1: Estimate Darkest hour
[totalIntensity,timeArr]=estimate_darkest_frame(fileStr);
[val,indx]=min(totalIntensity);
timeStr = datestr(timeArr(indx));
time = timeArr(indx);
dasc.sensorLoc = h5read(fileStr,'/DASC/sensorloc');

polarisIndx = find(stars.HIP==11767); % Validation star : Polaris

% Approximation
[stars.az,stars.el] = RADec2AzEl(rad2deg(stars.RA),rad2deg(stars.DEC),...
    dasc.sensorLoc(1),dasc.sensorLoc(2),datestr(time,'yyyy/mm/dd HH:MM:ss'));

% Sophesticated calculation of Polaris
[polaris.azObs,polaris.elObs] = get_star_az_el...
    (stars.RA(polarisIndx),stars.DEC(polarisIndx),...
    stars.pmRA(polarisIndx),stars.pmDEC(polarisIndx),stars.parallax(polarisIndx),...
    stars.RV(polarisIndx),time,deg2rad(dasc.sensorLoc(1)),deg2rad(dasc.sensorLoc(2)),dasc.sensorLoc(3));

%% Reading DASC data
asiPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC\20080326\';
asi1File = 'PKR_DASC_0000_20080326_103958.000.FITS';
ASI1 = fitsread([asiPath,asi1File]);
dasc.ASI = read_h5_variable_at_time_v2(fileStr,'/DASC/ASI',[],timeStr);
dasc.lat = (h5read(fileStr,'/DASC/lat'))';
dasc.lon = (h5read(fileStr,'/DASC/lon'))';

%% Project stars to 110 km 
projectionAltitude = 110-dasc.sensorLoc(3)./1000;
RE = 6.371*10^6;
stars.slantRange = -RE.*sind(stars.el) + modSign(stars.el).*sqrt((RE^2).*(sind(stars.el)).^2 + projectionAltitude.*1000.*(projectionAltitude.*1000+2.*RE)); 
polaris.slantRange = -RE.*sind(polaris.elObs) + modSign(polaris.elObs).*sqrt((RE^2).*(sind(polaris.elObs)).^2 + projectionAltitude.*1000.*(projectionAltitude.*1000+2.*RE)); 

starFilter = stars.el<10;
[stars.lat, stars.lon, stars.alt] = aer2geodetic(stars.az,stars.el,stars.slantRange,...
    dasc.sensorLoc(1),dasc.sensorLoc(2),dasc.sensorLoc(3),wgs84Ellipsoid('m'));
stars.lat(starFilter) = nan;
stars.lon(starFilter) = nan;

[polaris.lat, polaris.lon, polaris.alt] = aer2geodetic(polaris.azObs,polaris.elObs,polaris.slantRange,...
    dasc.sensorLoc(1),dasc.sensorLoc(2),dasc.sensorLoc(3),wgs84Ellipsoid('m'));

%% plot dasc given lat lon coordinates
hFig=figure;
p=create_panels(hFig,'totalPanelNo',2,'panelHeight',100,...
    'panelBreadth',100,'demargin',1);
p(1,1).select();
latLim=[63 67];
lonLim = [-153 -143];
deltaLat = 1;
deltaLon = 5;
ax1=axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);
axis off
load coastlines
plotm(coastlat,coastlon)
hold on;

h=pcolorm(latCal,lonCal,ASI1); 
set(h,'EdgeColor','none');
title([timeStr,' Toshi Calibration']);

hold on;
scatterm(stars.lat,stars.lon,10*stars.relIntensity,'r');
% plot_DASC_geodetic(ASI1,time,latCal,lonCal,1024,[60 70],[-160 -135],2,10,'test');

% Vega
starIndx = find(contains(deblank(stars.name),'Vega','IgnoreCase',true));
scatterm(stars.lat(starIndx),stars.lon(starIndx),10,'y');
textm(stars.lat(starIndx),stars.lon(starIndx),' Vega','Color','y');
scatterm(polaris.lat,polaris.lon,10,'c');
textm(polaris.lat,polaris.lon,' Polaris','Color','c');

colormap(get_colormap('k','w'));
colorbar('eastoutside');
caxis([300,500]);
%
p(1,2).select();
latLim=[63 67];
lonLim = [-153 -143];
deltaLat = 1;
deltaLon = 5;
ax2=axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);
axis off
load coastlines
plotm(coastlat,coastlon)
hold on;

h1=pcolorm(dasc.lat,dasc.lon,dasc.ASI); 
set(h1,'EdgeColor','none');
title([timeStr,' before 2011 calibration']);

hold on;
scatterm(stars.lat,stars.lon,10*stars.relIntensity,'r');
% plot_DASC_geodetic(ASI1,time,latCal,lonCal,1024,[60 70],[-160 -135],2,10,'test');

% Vega
starIndx = find(contains(deblank(stars.name),'Vega','IgnoreCase',true));
scatterm(stars.lat(starIndx),stars.lon(starIndx),10,'y');
textm(stars.lat(starIndx),stars.lon(starIndx),' Vega','Color','y');
scatterm(polaris.lat,polaris.lon,10,'c');
textm(polaris.lat,polaris.lon,' Polaris','Color','c');
colormap(get_colormap('k','w'));
colorbar('eastoutside');
caxis([300,500]);

%% Functions
function [totalIntensity,timeArr]=estimate_darkest_frame(h5FileStr)

asi = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
timeArr = unixtime2matlab((h5read(h5FileStr,'/DASC/time'))');
totalIntensity = sum(sum(asi,3),2);

end

function y = modSign(x)
    for i=1:1:length(x)
        if x(i)==0
            x(i) = +1;
        end
    end
    y = sign(x);
end
