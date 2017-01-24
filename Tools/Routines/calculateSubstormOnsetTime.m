%% Calculate the substorm onset time, and plot a figure
% And tag several other events
clear all;

load '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/thmasi_data_26_Mar_2008.mat'
load '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/pfisr_data_26_Mar_2008.mat'
load '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat'
load '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Global/GeomagneticIndex_data_26_Mar_2008.mat'

asi.fykn.intensitySum=sum(asi.fykn.intensity,2);
asi.gako.intensitySum=sum(asi.gako.intensity,2);
pfisrBeamAvg.density.electronSum = sum(pfisrBeamAvg.density.electron,1);

totalPanelNo=4;

hFig=figure(1);
resize_figure(hFig);

clf
p=panel();
p.pack(1);

panelSize = 50; %in mm
panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end;
p(1).pack(panelDefinition);
p.marginleft=30;
p.marginright=10;
p(1).de.margin=4;
% p.fontsize=12;
p.select('all');
timeMinStr='2008-03-26/08:05:00';
timeMaxStr='2008-03-26/12:55:00';
timeTick=0.5;


p(1,1).select();
plot(geomagIndex.time,geomagIndex.al,'k');
label_time_axis(geomagIndex.time, false ,0.5, timeMinStr, timeMaxStr);
ylabel({'AL Index', '[\muT]'});
grid on;

p(1,2).select();
plot(thd.magneticVectorFieldFGM.time,thd.magneticVectorFieldFGM.Bz,'k');
hold on;
plot(thd.magneticVectorFieldFGM.time,thd.magneticVectorFieldFGM.Bx,'--r');
label_time_axis(thd.magneticVectorFieldFGM.time, false ,0.5, timeMinStr, timeMaxStr);
ylabel({'thd B-field', '[nT]'});
legend('Bz', 'Bx','Orientation','horizontal','Location','southwest');
grid on;

p(1,3).select();
plot(asi.fykn.time,asi.fykn.intensitySum,'k');
hold on;
plot(asi.gako.time,asi.gako.intensitySum,'--r');
set(gca,'YScale','linear');
label_time_axis(asi.gako.time, false ,0.5, timeMinStr, timeMaxStr);
legend('FYKN', 'GAKO','Orientation','horizontal','Location','northwest');
ylabel({'Total Intensity', 'a.u.'});
grid on;

p(1,4).select();
plot(pfisrBeamAvg.density.time,pfisrBeamAvg.density.electronSum,'k');
label_time_axis(pfisrBeamAvg.density.time, true ,0.5, timeMinStr, timeMaxStr);
ylabel({'Total electron','density', '[m^-^3]'});
grid on;

