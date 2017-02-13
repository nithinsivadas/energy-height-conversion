clear all;
allData = get_all_time_series_data();

data = allData(4);

fileNameStr='/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
[datamag.zValue, datamag.yAxis, datamag.time] =...
    get_pfisr_variable(fileNameStr, 'popl', -1);
datamag.zValue=interp_nans(datamag.zValue);
%%
timeMin = datenum('26-Mar-2008 10:00');
timeMax = datenum('26-Mar-2008 11:00');

timeAvgDataAvg_1 = get_time_avg_time_series_data(data.zValue, data.time, timeMin, timeMax);
timeAvgDataMag_1 = get_time_avg_time_series_data(datamag.zValue, datamag.time, timeMin, timeMax);

timeMin = datenum('26-Mar-2008 11:00');
timeMax = datenum('26-Mar-2008 11:30');

timeAvgDataAvg_2 = get_time_avg_time_series_data(data.zValue, data.time, timeMin, timeMax);
timeAvgDataMag_2 = get_time_avg_time_series_data(datamag.zValue, datamag.time, timeMin, timeMax);

timeMin = datenum('26-Mar-2008 11:44');
timeMax = datenum('26-Mar-2008 11:48');

timeAvgDataAvg_3 = get_time_avg_time_series_data(data.zValue, data.time, timeMin, timeMax);
timeAvgDataMag_3 = get_time_avg_time_series_data(datamag.zValue, datamag.time, timeMin, timeMax);


%%
hFig=figure(5);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 100; %in mm
p.pack({{panelSize} {panelSize}},{{80}});
p.marginleft=30;
p.marginright=10;
p.fontsize=12;
p.de.margin=30;
p(1).margintop=60;
p(2).margintop=30;
p.select('all');

p(1,1).select();
plot(timeAvgDataAvg_1,data.yAxis,'-k'); 
hold on;
plot(timeAvgDataAvg_2,data.yAxis,'-r'); 
hold on;
plot(timeAvgDataAvg_3,data.yAxis,'-g'); 
set(gca, 'XScale','log');
legend('10:00-11:00 UT','11:00-11:30 UT','11:44-11:48','Location','northwest');
ylabel('Altitude [km]');
xlabel({'Electron density [m^-^3]','Averaged along all beams'});
grid on;
ylim([50 200]);
set(gca, 'YTick',[50, 60, 70, 80, 90, 100, 110, 120, 150, 200]);

p(2,1).select();
plot(timeAvgDataMag_1,datamag.yAxis,'-k'); 
hold on;
plot(timeAvgDataMag_2,datamag.yAxis,'-r'); 
hold on;
plot(timeAvgDataMag_3,datamag.yAxis,'-g'); 
set(gca, 'XScale','log');
legend('10:00-11:00 UT','11:00-11:30 UT','11:44-11:48','Location','northwest');
ylabel('Altitude [km]');
xlabel({'Electron density [m^-^3]','Along magnetic-field aligned beam'});
grid on;
ylim([50 200]);
set(gca, 'YTick',[50, 60, 70, 80, 90, 100, 110, 120, 150, 200]);