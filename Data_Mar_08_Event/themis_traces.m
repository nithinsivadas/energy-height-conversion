% Create THEMIS Traces with  PFISR ground location
clear all;

% Creating Coastlines
figure('Color','w') 
axesm('lambertstd','MapLatLimit',[40 80],'MapLonLimit',[-200 -120])
axis off; framem on; gridm on; mlabel on; plabel on;
load coastlines
plotm(coastlat,coastlon) 

load('thd_ntrace.mat');

%Generating a time array
for i=1:1:length(thd_ntrace.time)
    time(i) = datenum(strcat('26-Mar-2008/',thd_ntrace.time(i)));
end;

% Decinding time indices that require a tick mark
index(1)=1;
for i=2:1:ceil((max(time)-min(time))/(1/24))
[c index(i)] = min(abs(time - (datenum('26-Mar-2008/08:00:00')+(1/24)*(i-1))));
end;
% index(i+1)=length(time);

% Plotting THEMIS Trace
plotm(thd_ntrace.lat,thd_ntrace.long,'g-');
textm(thd_ntrace.lat(index)-2, thd_ntrace.long(index)-4, datestr(time(index),'HH:MM'),'FontSize',12); 
plotm(thd_ntrace.lat(index), thd_ntrace.long(index), '+black'); 
textm(thd_ntrace.lat(index(end))+2, thd_ntrace.long(index(end))-8, 'THEMIS-D Track','Color','green','FontSize',14); 

pfisrlat=65.12;
pfisrlong=-147.48;
plotm(pfisrlat,pfisrlong,'r*');
textm(pfisrlat+1, pfisrlong, ' Poker Flat','Color','red','FontSize',14) 
% datestr(datenum(strcat('26-Mar-2008/',thd_ntrace.time(1))),'HH:MM')]]

lat1 = thd_ntrace.lat*pi/180;
lat2 = pfisrlat*pi/180;
lon1 = thd_ntrace.long*pi/180;
lon2 = pfisrlong*pi/180;
dlat = lat2-lat1;
dlon = lon2-lon1;

a = (sin(dlat/2).^2)+cos(lat1).*cos(lat2).*(sin(dlon/2)).^2;
c = 2.*atan2(a.^0.5, (1-a).^0.5);
d = 6371.*c;

[M,I]=min(d);
display(['Minimum distance between PFISR and THM-D is ',num2str(M),' km achieved at ',datestr(time(I)),' UT']);

