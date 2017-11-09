function [ figureHandle ] = plot_2D_integrated_conductivity_map...
    (data, mode, timeNo, setMapOn, latLim, lonLim)
%plot_2D_integrated_conductivity_map Plotting a 2D map of integrated conductivity
% Input
% data : data.sigma_P --> Matrix of electron density point at each[position, time]
%      : the measurements can be organized in any order
%      : data.time --> Time index number
%      : data.lat -->
%      : data.lon
%      : data.projectionAltitude--> altitude slice of the map [km]
%      : data.time --> matlab datenum time of the map
% setMapOn: 'true' value projects on a map
% latLim : the min and max value of the latitude to plot
% lonLim : the min and max value of the longitude to plot

%% Setting default values

if nargin <4
    setMapOn=true;
end

%% Initializing coordinates
lat = data.lat;
lon = data.lon;

%% Generating data slice

% Loading default values for the lat and lon limit, if not specified
if nargin<5
latLim = [min(lat(:)) max(lat(:))];
end
if nargin<6
lonLim = [min(lon(:)) max(lon(:))];
end
%% Plotting map
imageSize=512;
latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
if mode=='H'
    v=data.sigma_H(:,timeNo);
elseif mode=='P'
    v=data.sigma_P(:,timeNo);
elseif mode=='R'
    v=data.ratio(:,timeNo);
end
[Xlatq,Ylonq] = meshgrid(latq,lonq);
Vq = griddata(lat,lon,v,Xlatq(:),Ylonq(:),'nearest');
Vq1= reshape(Vq,512,512);

% Ensuring that points outside the data field are not plotted
K=convhull(lat,lon);
in = inpolygon(Xlatq,Ylonq,lat(K),lon(K));
Vq2 = Vq1.*in;
Vq2(Vq2==0)=nan;
if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1)) (latLim(2))],...
        'MapLonLimit',[(lonLim(1)) (lonLim(2))],'FontSize',6,...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',0.5,'MLineLocation',1,'MLabelRound',-1,'PLabelRound',-1);
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(Xlatq.*in,Ylonq.*in,(Vq2));
    hold on;
else
    h2=pcolor(lonq,latq,log10(Vq));
%     hold on;
%     title(['PFISR: ', num2str(altitudeSelected),' km ',datestr(data.time(timeNo),'HH:MM:SS'),' UT']);
end

set(h2,'EdgeColor','none');

end
