%% Testing a star chart

stars=load('stars');
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
%%
timeStr = '26 Mar 2008 10:30:50';
time = datenum(timeStr);
inCoord = [stars.long,stars.lat];
% sensorLoc = [+65.126+0.1,-147.4789+1,0.689];
sensorLoc = [62.407-0.1,-145.1580+0.9,0];
% outCoord = horiz_coo(inCoord,jd,deg2rad(sensorLoc),'h');
[stars.Az,stars.El] = RADec2AzEl(stars.long,stars.lat,datestr(time,'yyyy/mm/dd HH:MM:ss'),sensorLoc(1),sensorLoc(2),sensorLoc(3));
stars.El = stars.El;
stars.El(stars.El<1)=nan;

projectionAltitude = 90; % km
RE = 6.371*10^6;
projectionAltitude1 = (projectionAltitude - sensorLoc(3))*1000;
stars.slantRange = projectionAltitude1./sind(stars.El);
% stars.slantRange = -RE.*sind(stars.El) + modSign(stars.El).*sqrt((RE^2).*(sind(stars.El)).^2 + projectionAltitude1.*(projectionAltitude1+2.*RE)); 

[stars.glatitude,stars.glongitude,stars.galtitude] = aer2geodetic(stars.Az,stars.El,stars.slantRange,sensorLoc(1),sensorLoc(2),sensorLoc(3)*1000,wgs84Ellipsoid('m'));

% latLim = [58,69];
% lonLim = [-160,-130];
% deltaLat = 2;
% deltaLon = 10;

starMag = 2;
starFilter = stars.vmag<=starMag;
latLim = [60,65];
lonLim = [-153,-137];
deltaLat = 2;
deltaLon = 10;

storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\Figure1b\';
h=figure('visible','on');

[ax1]=combine_2D_plots_v2(fileStr,h,...
    'maps',{'OpticalImage'},...
    'sites',{'gako'},...
    'thisTime',time(1),...
    'latLim',latLim,...
    'lonLim',lonLim,...
    'elCutOff',20,...
    'deltaLat',deltaLat,...
    'deltaLon',deltaLon,...
    'opticalLim',[250 450],...
    'peakIonizationAltitude',85,...
    'setStoreImage',false);
cm=get_colormap('w',[0,0.2,0]);
colormap(cm);
hold on;
scatterm(stars.glatitude(starFilter),stars.glongitude(starFilter),...
    stars.relintensity(starFilter)*100,'r');

function y = modSign(x)
    for i=1:1:length(x)
        if x(i)==0
            x(i) = +1;
        end
    end
    y = sign(x);
end