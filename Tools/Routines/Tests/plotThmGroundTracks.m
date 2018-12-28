%% Satellite tracks on ground
%%
% load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS89\lossConeFluxTHMADE_20080326_TS89_with_footpoints.mat');
%%
thisTimeStr = '26 Mar 2008 11:10';
omnih5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
%%
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
h=figure;
combine_2D_plots_v2(fileStr,h,...
    'maps',{'OpticalImage','OpticalImage','OpticalImage','OpticalImage','OpticalImage','OpticalImage','EnergyFluxMap'},...
    'sites',{'fykn','gako','whit','mcgr','kian','inuv','pokerFlat'},...
    'thisTime',datenum(thisTimeStr),...
    'latLim',[58,72],...
    'lonLim',[-179,-120],...
    'elCutOff',15,...
    'deltaLat',2,...
    'deltaLon',10,...
    'opticalLim',[400 550],...
    'peakIonizationAltitude',85,...
    'setStoreImage',false);

% combine_2D_plots_v2(fileStr,h,...
%     'maps',{'OpticalImage','EnergyFluxMap'},...
%     'sites',{'pokerFlat','pokerFlat'},...
%     'thisTime',datenum(thisTimeStr),...
%     'opticalLim',[400 550],...
%     'elCutOff',15,...
%     'latLim',[63,67],...
%     'lonLim',[-152 -142],...
%     'peakIonizationAltitude',85,...
%     'setStoreImage',false);
%% Calculating THM Foot prints
hold on;
timeMinStr = '26 Mar 2008 08:00';
timeMaxStr = '26 Mar 2008 13:00';
% kext = find_irbem_magFieldModelNo('TS96');
% options = [0,0,0,0,0];
% sysaxes = 2;
% stop_alt = 110;
% [maginput,timeMaginput] = generate_maginput(omnih5,timeMinStr,timeMaxStr);
% maginput=interp_nans(maginput);
% maginput = filter_irbem_maginput(kext,maginput);
% thmData=process_themis_data('26 Mar 2008');
%%
% tic
% 
% tRange = find_time(thmData.thd.state.time,timeMinStr):1:find_time(thmData.thd.state.time,timeMaxStr);
% XFoot = zeros(length(tRange),3);
% XFootG = zeros(length(tRange),3);
% for i = 1:1:length(tRange)
%     thisTime = thmData.thd.state.time(tRange(i));
%     thisMaginput = interp1(timeMaginput,maginput,thisTime,'nearest');
%     x1=thmData.thd.state.XYZ_GSM(tRange(i),1);
%     x2=thmData.thd.state.XYZ_GSM(tRange(i),2);
%     x3=thmData.thd.state.XYZ_GSM(tRange(i),3);
%     XFoot(i,:)=onera_desp_lib_find_foot_point(kext,...
%         options,sysaxes,thisTime,x1,x2,x3,stop_alt,+1,thisMaginput);
%         XFoot(i,:)=onera_desp_lib_find_foot_point(kext,...
%         options,sysaxes,thisTime,x1,x2,x3,stop_alt,+1,thisMaginput);
%     XFootG(i,:)=geopack_find_foot_point(kext,[],...
%             2,thisTime,...
%             x1,...
%             x2,...
%             x3,...
%             stop_alt,+1,thisMaginput);
% end
% padData.thd.time = thmData.thd.state.time(tRange);
% padData.thd.latFoot = XFootG(:,2);
% padData.thd.lonFoot = XFootG(:,3);
% toc
% 
% tic
% tRange = find_time(thmData.the.state.time,timeMinStr):1:find_time(thmData.the.state.time,timeMaxStr);
% XFoot = zeros(length(tRange),3);
% XFootG = zeros(length(tRange),3);
% for i = 1:1:length(tRange)
%     thisTime = thmData.the.state.time(tRange(i));
%     thisMaginput = interp1(timeMaginput,maginput,thisTime,'nearest');
%     x1=thmData.the.state.XYZ_GSM(tRange(i),1);
%     x2=thmData.the.state.XYZ_GSM(tRange(i),2);
%     x3=thmData.the.state.XYZ_GSM(tRange(i),3);
%     XFoot(i,:)=onera_desp_lib_find_foot_point(kext,...
%         options,sysaxes,thisTime,x1,x2,x3,stop_alt,+1,thisMaginput);
%     XFootG(i,:)=geopack_find_foot_point(kext,[],...
%             2,thisTime,...
%             x1,...
%             x2,...
%             x3,...
%             stop_alt,+1,thisMaginput);
% end
% 
% padData.the.time = thmData.the.state.time(tRange);
% padData.the.latFoot = XFootG(:,2);
% padData.the.lonFoot = XFootG(:,3);
% toc
%%
% lineWidth = 1;
% 
% % pfisr
% hold on;
% m=plotm(65.12,-147.47,'or');
% % set(m, {'MarkerFaceColor'}, 'r'); 
% % thd
% time = padData.thd.time;
% lat = conv(padData.thd.latFoot,ones(5,1)/5,'same');
% lon = conv(padData.thd.lonFoot,ones(5,1)/5,'same');
% 
% timeMinIndx = find_time(time,timeMinStr);
% timeMaxIndx = find_time(time,timeMaxStr);
% thisTimeIndx = find_time(time,thisTimeStr);
% tIndx = timeMinIndx:timeMaxIndx;
% 
% plotm(lat(tIndx),lon(tIndx),'-b','LineWidth',lineWidth);
% plotm(lat(thisTimeIndx),lon(thisTimeIndx),'*k');
% hhmm = {'09:30','09:45','10:00','10:15','10:30','10:45','11:00'};
% tArray = datenum(strcat('26 Mar 2008'," ",hhmm));
% plot_time_markers(time,lat,lon,tArray,'Southwest','k',[0.8 0.8 0.8]);
% 
% % the
% time = padData.the.time;
% lat = conv(padData.the.latFoot,ones(5,1)/5,'same');
% lon = conv(padData.the.lonFoot,ones(5,1)/5,'same');
% 
% timeMinIndx = find_time(time,timeMinStr);
% timeMaxIndx = find_time(time,timeMaxStr);
% thisTimeIndx = find_time(time,thisTimeStr);
% tIndx = timeMinIndx:timeMaxIndx;
% 
% plotm(lat(tIndx),lon(tIndx),'-c','LineWidth',lineWidth);
% plotm(lat(thisTimeIndx),lon(thisTimeIndx),'*k');
% % hhmm = {'09:00','10:00','11:00'};
% hhmm = {'08:45','09:00','09:15','09:30'};
% tArray = datenum(strcat('26 Mar 2008'," ",hhmm));
% plot_time_markers(time,lat,lon,tArray,'Northeast','k',[0.8 0.8 0.8]);
% 
% % tha
% time = padData.tha.time;
% lat = conv(padData.tha.latFoot,ones(5,1)/5,'same');
% lon = conv(padData.tha.lonFoot,ones(5,1)/5,'same');
% 
% timeMinIndx = find_time(time,timeMinStr);
% timeMaxIndx = find_time(time,timeMaxStr);
% thisTimeIndx = find_time(time,thisTimeStr);
% tIndx = timeMinIndx:timeMaxIndx;
% 
% plotm(lat(tIndx),lon(tIndx),'-g','LineWidth',lineWidth);
% plotm(lat(thisTimeIndx),lon(thisTimeIndx),'*k');
% hhmm = {'09:00','10:00','11:00','11:30'};
% tArray = datenum(strcat('26 Mar 2008'," ",hhmm));
% plot_time_markers(time,lat,lon,tArray,'Northeast','k',[0.8 0.8 0.8]);

%%


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