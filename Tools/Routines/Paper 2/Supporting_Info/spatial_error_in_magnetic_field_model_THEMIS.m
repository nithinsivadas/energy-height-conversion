%% THEMIS-D mapping error
%% Load Data
load('G:\My Drive\Research\Projects\Collaborations\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
%%
omnih5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
% fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
% poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
% load(poesMatFile);
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');

%%
% magFieldNo = find_irbem_magFieldModelNo('TS01');
% maginput = filter_irbem_maginput(magFieldNo,maginput);
thd = padData.thd;
%%
% POES trajectory
timeMinStr = '26 Mar 2008 08:00:00';
timeMaxStr = '26 Mar 2008 12:00:00';
timeIndx = find_time(thd.time,timeMinStr):1:find_time(thd.time,timeMaxStr);
thdtime = thd.time(timeIndx);
thdGSE = thd.XYZ_GSE(timeIndx,:);

kext1 = 1:1:11;
for k=1:1:length(kext1)
    maginput1 = filter_irbem_maginput(kext1(k),maginput);
    for i = 1:1:length(thdtime)
    thisMaginput = interp1(timeMaginput,maginput1,thdtime(i));
    [~,thdEqPoint1(i,k,:)] = onera_desp_lib_find_magequator(kext1(k),[0,0,0,0,0],3,thdtime(i),...
        thdGSE(i,1),thdGSE(i,2),thdGSE(i,3),thisMaginput);
    thdEqGSM(i,k,:) = onera_desp_lib_coord_trans(thdEqPoint1(i,k,:),...
        [1 2],thdtime(i));
%     thdNFoot1(i,k,:)=onera_desp_lib_find_foot_point(kext1(k),[0,0,0,0,0],...
%         1,thdtime(i),thdEqPoint1(i,1),thdEqPoint1(i,2),thdEqPoint1(i,3),85,+1,thisMaginput);
    thdNFoot1(i,k,:)=onera_desp_lib_find_foot_point(kext1(k),[0,0,0,0,0],...
        3,thdtime(i),thdGSE(i,1),thdGSE(i,2),thdGSE(i,3),85,+1,thisMaginput);
    end
end

for iTime = 1:1:length(thdtime)
    sigmaLatFP(iTime) = nanstd(thdNFoot1(iTime,:,2));
    meanLatFP(iTime) = nanmean(thdNFoot1(iTime,:,2));
    nFP(iTime) = sum(~isnan(thdNFoot1(iTime,:,2)));
    sigmaLonFP(iTime) = nanstd(thdNFoot1(iTime,:,3));
    meanLonFP(iTime) = nanmean(thdNFoot1(iTime,:,3));
    
    sigmaMagEq(iTime,:) = [nanstd(thdEqGSM(iTime,:,1)),nanstd(thdEqGSM(iTime,:,2)),nanstd(thdEqGSM(iTime,:,3))];
    meanMagEq(iTime,:) = [nanmean(thdEqGSM(iTime,:,1)),nanmean(thdEqGSM(iTime,:,2)),nanmean(thdEqGSM(iTime,:,3))];
    nEq(iTime) = sum(~isnan(thdEqGSM(iTime,:,1).*thdEqGSM(iTime,:,2).*thdEqGSM(iTime,:,3)));
    
end


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
end
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
boxplot(thdEqGSM(:,:,1)','positions',thdtime,'Symbol','');
SEMagEqx = 1.96*sigmaMagEq(:,1)./sqrt(nEq)';
% plot(thdtime,meanMagEq(:,1),'-k');
% hold on;
% plot(thdtime,meanMagEq(:,1)+1.96*sigmaMagEq(:,1)./sqrt(nEq)','-r');
% plot(thdtime,meanMagEq(:,1)-1.96*sigmaMagEq(:,1)./sqrt(nEq)','-b');

label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Southwest');
ylabel({'Themis-D','Mag. Equator','X_G_S_M [R_E]'});
ylim([-13 -8]);

q(2).select();

boxplot(thdEqGSM(:,:,3)','positions',thdtime,'Symbol','');
SEMagEqz = 1.96*sigmaMagEq(:,3)./sqrt(nEq)';
% plot(thdtime,meanMagEq(:,3),'-k');
% hold on;
% plot(thdtime,meanMagEq(:,3)+1.96*sigmaMagEq(:,3)./sqrt(nEq)','-r');
% plot(thdtime,meanMagEq(:,3)-1.96*sigmaMagEq(:,3)./sqrt(nEq)','-b');

label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Northwest');
ylabel({'Themis-D','Mag. Equator','Z_G_S_M [R_E]'});
ylim([-0.5 1]);

q(3).select();
boxplot(thdNFoot1(:,:,2)','positions',thdtime,'Symbol','');
SELatFP = 1.96*sigmaLatFP./sqrt(nFP);
% plot(thdtime,meanLatFP,'-k');
% plot(thdtime,meanLatFP+1.96*sigmaLatFP./sqrt(nFP),'-r');
% plot(thdtime,meanLatFP-1.96*sigmaLatFP./sqrt(nFP),'-b');

label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Northwest');
ylabel({'Themis-D','North FP','Lat[°N]'});

q(4).select();
boxplot(thdNFoot1(:,:,3)','positions',thdtime,'Symbol','');
SELonFP = 1.96*sigmaLonFP./sqrt(nFP);
% plot(thdtime,meanLonFP,'-k');
% hold on;
% plot(thdtime,meanLonFP+1.96*sigmaLonFP./sqrt(nFP),'-r');
% plot(thdtime,meanLonFP-1.96*sigmaLonFP./sqrt(nFP),'-b');

label_time_axis(thdtime,true,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Southwest');
ylabel({'Themis-D','North FP','Lon[°E]'});

%%
totalPanelNo=4;

hFig=figure(2);

clf
p=panel();
p.pack(1);

panelSize = 35; %in mm
demargin = 4;
panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end
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
plot(thdtime,SEMagEqx);
grid on;
label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
ylabel({'Themis-D','Mag. Equator','X_G_S_M S.E. [R_E]'});
% ylim([-13 -8]);

q(2).select();

plot(thdtime,SEMagEqz);
grid on;
% plot(thdtime,meanMagEq(:,3),'-k');
% hold on;
% plot(thdtime,meanMagEq(:,3)+1.96*sigmaMagEq(:,3)./sqrt(nEq)','-r');
% plot(thdtime,meanMagEq(:,3)-1.96*sigmaMagEq(:,3)./sqrt(nEq)','-b');

label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Northwest');
ylabel({'Themis-D','Mag. Equator','Z_G_S_M S.E. [R_E]'});
% ylim([-0.5 1]);

q(3).select();
plot(thdtime,SELatFP);
grid on;
% plot(thdtime,meanLatFP,'-k');
% plot(thdtime,meanLatFP+1.96*sigmaLatFP./sqrt(nFP),'-r');
% plot(thdtime,meanLatFP-1.96*sigmaLatFP./sqrt(nFP),'-b');

label_time_axis(thdtime,false,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Northwest');
ylabel({'Themis-D','North FP','S.E. Lat[°N]'});

q(4).select();
plot(thdtime,SELonFP);
% SELonFP = 1.96*sigmaLonFP./sqrt(nFP);
% plot(thdtime,meanLonFP,'-k');
% hold on;
% plot(thdtime,meanLonFP+1.96*sigmaLonFP./sqrt(nFP),'-r');
% plot(thdtime,meanLonFP-1.96*sigmaLonFP./sqrt(nFP),'-b');
grid on;
label_time_axis(thdtime,true,timeTick,timeMinStr,timeMaxStr);
% legend('Mean','Upper 95% limit','Lower 95% limit','Location','Southwest');
ylabel({'Themis-D','North FP','S.E. Lon[°E]'});



function plot_time_markers(time,lat,lon,markTimeArr,orientation,tickColor,textColor)
    if nargin<7
        textColor = 'k';
    end
    if nargin <6
        tickColor = 'k';
    end

    if strcmpi(orientation(1:5),'North')
%         dlat = +0.2;
        dlat = +0.05;
    elseif strcmpi(orientation(1:5),'South')
%         dlat = -0.2;
        dlat = -0.05;
    end
    if strcmpi(orientation(6:9),'East')
%         dlon = +0.5;
        dlon = +0.1;
    elseif strcmpi(orientation(6:9),'West')
%         dlon = -3;
        dlon = -1;
    end

    for i = 1:length(markTimeArr)
        thisTimeIndx = find_time(time,datestr(markTimeArr(i)));
        plotm(lat(thisTimeIndx),lon(thisTimeIndx),'.','Color',tickColor);
        textm(lat(thisTimeIndx)+dlat,lon(thisTimeIndx)+dlon...
            ,datestr(markTimeArr(i),'HH:MM'),...
            'Color',textColor);
                %         create_perp_line(lat,lon,thisTimeIndx,orientation(1:5));
    end

end

function create_perp_line(lat,lon,markTimeArrIndx,orientation,tickLength)
    if nargin<6
        tickLength = 0.1;
    end

    if strcmpi(orientation(1:5),'North')
        dlat = 1;
    elseif strcmpi(orientation(1:5),'South')
        dlat = -1;
    end

    for i = length(markTimeArrIndx)
        thisTimeIndx = markTimeArrIndx(i);
        lat1 = lat(thisTimeIndx);
        lon1 = lon(thisTimeIndx);
        dxlon = mean(diff(lon(thisTimeIndx-1:thisTimeIndx+1)));
        dylat = gradient(lat(thisTimeIndx-1:thisTimeIndx+1),dxlon);
        slope = -1./dylat(2);
        c = lat1-slope*lon(1);
        lat2 = lat1+dlat*tickLength;
        lon2 = (lat2-c)./slope;
        hold on;
        plotm([lat1,lat2],[lon1,lon2],'k');
    end

end

function [fov] = create_dasc_fov(fileStr, elCutOff, projectionAltitude)

    dascData = get_2D_plot_inputs_time_independent(fileStr,'plotModeStr','OpticalImage','site','pokerFlat');
    az = 0:0.1:360;
    el = elCutOff*ones(size(az));
    alt = projectionAltitude;
    C = define_universal_constants;
    slantRange = (-C.RE.*sind(el) + sqrt((C.RE^2).*(sind(el)).^2 + alt.*1000.*(alt.*1000+2.*C.RE)))./1000; 
    
%     slantRange = alt./sind(el);
    [fov.lat,fov.lon,fov.alt] = aer2geodetic(az,el,slantRange,dascData.sensorloc(1),dascData.sensorloc(2),dascData.sensorloc(3)./1000,wgs84Ellipsoid('km'));
    
end