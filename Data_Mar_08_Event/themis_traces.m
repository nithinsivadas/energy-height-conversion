% Create THEMIS Traces with  PFISR ground location

figure('Color','w') 
axesm('lambertstd','MapLatLimit',[40 80],'MapLonLimit',[-200 -120])
axis off; framem on; gridm on; mlabel on; plabel on;
load coastlines
plotm(coastlat,coastlon) 

pfisrlat=65;
pfisrlong=-147.5;
plotm(pfisrlat,pfisrlong,'r*');
textm(pfisrlat, pfisrlong, ' Poker Flat') 

load('thd_ntrace.mat');
plotm(thd_ntrace.lat,thd_ntrace.long,'g-');
textm(thd_ntrace.lat([1,293]), thd_ntrace.long([1,293]), thd_ntrace.time([1,293])) 



