%% Themis Test
% created by J. Brent Parham 4/17/2018
clc; clear all;
load auroraColor
for i=0:120
figure(1); clf;
    worldmap([50 70], [-135 -100]);
    geoshow('landareas.shp','FaceColor',[.9 .9 .9])
    hold on;
    
    dateThem=datetime('14-Jan-2018 4:40:00') + seconds(5*i);
    
    camera=themisRead(dateThem);
    pcolorm(camera(1).lat,camera(1).lon,camera(1).themImage,'facealpha',0.8);
    pcolorm(camera(2).lat,camera(2).lon,camera(2).themImage,'facealpha',0.8);

    caxis([0 .15])
    title(datestr(dateThem))

    
    colormap(auroraColor)
    saveas(gcf,[sprintf('%03d',i),'storm.jpg'])
end