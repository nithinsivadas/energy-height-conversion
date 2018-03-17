timestr = ['21 Aug 2017 09:00';'21 Aug 2017 10:00';...
    '21 Aug 2017 11:00';...
    '21 Aug 2017 12:00';...
    '21 Aug 2017 13:00';
    '21 Aug 2017 14:00';
    '21 Aug 2017 15:00']; 
time = datenum(timestr);
lat = 40;
lon = -100;
altitudeIRI = 60:5:400;
coordinateSystem = 'geodetic';
curlDir='C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\External Tools';
figure;
for i=1:1:length(time)
    IRI=iri2012(time(i),lat,lon,altitudeIRI, true, coordinateSystem,curlDir, [], [], [], [], [], [], [], [], [], [] ,[] ,[], [], [], 'RBV10/TTS03' );
    Ni=IRI(:,1)';
    plot(Ni, altitudeIRI);
    hold on;
end
set(gca,'Xscale','log');
xlim([10^6 10^12]);
xlabel('Ion density [m^-^3]');
ylabel('Altitude [km]');
legend(timestr,'location','northwest');
title(['Lat = ',num2str(lat),' ^0',' Lon = ',num2str(lon),' ^0']);
