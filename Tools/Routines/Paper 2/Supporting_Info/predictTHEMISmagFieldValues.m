%% Predicting THEMIS-D magnetic field values using magnetic field models
% and comparing it with measured values. 
clear all;
%% Magnetic field model paramters
dateStr = '26 Mar 2008';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
magFieldModelNo = 7;
timeMinStr = '26 Mar 2008 08:00';
timeMaxStr = '26 Mar 2008 12:30';

%% Downloading themis state data
[themisData, matFilePath] = process_themis_data(dateStr);

%% Retrieving maginput
[maginput,timeMaginput]=generate_maginput(omniH5FileStr,...
    themisData.thd.state.time(1),themisData.thd.state.time(end));

%% Calculating magnetic field at themis-D location
magFieldModelStr=find_irbem_magFieldModelStr(magFieldModelNo);
maginput = filter_irbem_maginput(magFieldModelNo,maginput);
x1 = themisData.thd.state.XYZ_GEO(:,1);
x2 = themisData.thd.state.XYZ_GEO(:,2);
x3 = themisData.thd.state.XYZ_GEO(:,3);

timeMinIndx = find_time(timeMaginput,timeMinStr);
timeMaxIndx = find_time(timeMaginput,timeMaxStr);
timeArray = timeMinIndx: 1: timeMaxIndx;
multiWaitbar('Calculating',0);
dtime = 1./length(timeArray);
for itime = 1:length(timeArray)
    [Bgeo(itime,:),B(itime)] = onera_desp_lib_get_field(magFieldModelNo,[0 0 0 0],1,timeMaginput(1,timeArray(itime)),...
        x1(timeArray(itime)),x2(timeArray(itime)),x3(timeArray(itime)),maginput(timeArray(itime),:)); 
    multiWaitbar('Calculating','Increment',dtime);
end
multiWaitbar('Calculating','Close');

%% Download actual themis Bfield
themisDBfieldPath = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\THD_L2_FGM_6864.txt';
[themisTime,Bthm,Bxthm,Bythm,Bzthm] = get_themis_Bfield(themisDBfieldPath);
timeMinIndxThm = find_time(themisTime,timeMinStr);
timeMaxIndxThm = find_time(themisTime,timeMaxStr);
timeArrayThm = timeMinIndxThm:1:timeMaxIndxThm;

%%
totalPanelNo=4;

hFig=figure(1);

clf
p=panel();
p.pack(1);

panelSize = 35; %in mm
demargin = 4;
panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end;
p(1).pack(panelDefinition);

p.marginleft=35;
p.marginright=25;
p(1).de.margin=demargin;
% p.fontsize=12;
p.select('all');
timeMinStr='2008-03-26/08:05:00';
timeMaxStr='2008-03-26/12:00:00';
timeTick=0.5;
  
resize_figure(hFig, panelSize*totalPanelNo+4*(totalPanelNo+1)+20);

q=p(1);

q(1).select();
plot(timeMaginput(timeArray),Bgeo(:,1),'--k');
hold on;
plot(themisTime(timeArrayThm),Bxthm(timeArrayThm),'-r');
label_time_axis(timeMaginput(timeArray),false,timeTick,timeMinStr,timeMaxStr);
legend([magFieldModelStr,' Bx'],'THM-D Bx','Location','Southwest');
set(gca,'YLim',[-100 0]);
ylabel('[nT]');

q(2).select();
plot(timeMaginput(timeArray),Bgeo(:,2),'--k');
hold on;
plot(themisTime(timeArrayThm),Bythm(timeArrayThm),'-r');
label_time_axis(timeMaginput(timeArray),false,timeTick,timeMinStr,timeMaxStr);
legend([magFieldModelStr,' By'],'THM-D By','Location','Northwest');
set(gca,'YLim',[-50 50]);
ylabel('[nT]');

q(3).select();
plot(timeMaginput(timeArray),Bgeo(:,3),'--k');
hold on;
plot(themisTime(timeArrayThm),Bzthm(timeArrayThm),'-r');
label_time_axis(timeMaginput(timeArray),false,timeTick,timeMinStr,timeMaxStr);
legend([magFieldModelStr,' Bz'],'THM-D Bz','Location','Northwest');
set(gca,'YLim',[0 100]);
ylabel('[nT]');

q(4).select();
plot(timeMaginput(timeArray),B(:),'--k');
hold on;
plot(themisTime(timeArrayThm),Bthm(timeArrayThm),'-r');
label_time_axis(timeMaginput(timeArray),true,timeTick,timeMinStr,timeMaxStr);
legend([magFieldModelStr,' Bmag'],'THM-D Bmag','Location','Northwest');
set(gca,'YLim',[0 100]);
ylabel('[nT]');

%%
function [themisTime,B,Bx,By,Bz] = get_themis_Bfield(themisDBfieldPath)

    fileID = fopen(themisDBfieldPath); 
    formatSpec ='%s %s %f %f %f %f %f'; 
    C = textscan(fileID,formatSpec,'HeaderLines',64); 
    fclose(fileID);
    nLines=length(C{1});
    themisTime=datenum(join([C{1}(1:nLines-3),C{2}(1:nLines-3)]),'dd-mm-yyyy HH:MM:SS.FFF');
    B = C{3}(1:1:nLines-3);
    Bx = C{5}(1:1:nLines-3);
    By = C{6}(1:1:nLines-3);
    Bz = C{7}(1:1:nLines-3);
    
end   
