%% Saturation of AE given E_M
tic
dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
ace1 = extract_set_of_data(dataFolder,'ace');
wind1 = extract_set_of_data(dataFolder,'wind');
geotail1 = extract_set_of_data(dataFolder,'geotail');
omni = extract_omni_data(omniFile);

toc

%%
[~,~,wi]=intersect(geotail1.datetime,wind1.datetime);
[~,ai,wii]=intersect(ace1.datetime,wind1.datetime(wi));
[~,gi,~]=intersect(geotail1.datetime,wind1.datetime(wi(wii)));
wi = wi(wii);
[~,ai1,wi1]=intersect(ace1.datetime,wind1.datetime);
[~,gi1,wi2] = intersect(geotail1.datetime,wind1.datetime);
%%
ace = ace1(ai,:);
wind = wind1(wi,:);
geotail = geotail1(gi,:);

%%
ace_filter = ((ace.noseYGSE-ace.YGSE).^2 + (ace.noseZGSE-ace.ZGSE).^2).^0.5 < 10;
wind_filter = ((wind.noseYGSE-wind.YGSE).^2 + (wind.noseZGSE-wind.ZGSE).^2).^0.5 < 10;
geo_filter = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 10;

%%
figure;
minTime = [datetime(1999,08,06), datetime(1999,08,09)]; 
subplot(3,1,1)
plot(wind.datetime, wind.E_kl);
xlim(minTime);
subplot(3,1,2)
plot(geotail.datetime, geotail.E_kl);
xlim(minTime);
subplot(3,1,3)
plot(ace.datetime, ace.E_kl);
xlim(minTime);
%%
figure; 
fil=~isnan(geotail.E_kl) &  ~isnan(ace.E_kl) & ~isnan(wind.E_kl);
scatter(abs(wind.E_kl(fil)-geotail.E_kl(fil))./wind.E_kl(fil), abs(wind.E_kl(fil)-ace.E_kl(fil))./wind.E_kl(fil),10,'filled','MarkerFaceAlpha',0.2);
xlim([0 10]);
ylim([0 10]);

%% 
% Reponse of SML to E_m from OMNI database
filter = sqrt((wind.YGSE-wind.noseYGSE).^2 + (wind.ZGSE-wind.noseZGSE).^2)<10; 
Y = omni.Fsml(datenum(wind.datetime));
X = wind.E_kl.*10^-3;
Y(Y==999999)=nan;
[xindx, E] = discretize(X,100);

for i = 1:max(xindx)
    YgX(i) = nanmean(Y(xindx==i));
    NSamples(i) = sum(xindx==i & ~isnan(Y));
    SEM(i) = nanstd(Y(xindx==i))./sqrt(NSamples(i));
    ts(i,:) = tinv([0.025 0.975],NSamples(i)-1);
    CI(i,:) = YgX(i) + ts(i,:)*SEM(i);
    XBins(i) = 0.5*(E(i)+E(i+1));
end

%%

% Reponse of SML to E_m, when ACE and WIND measure almost the same thing
% filter = sqrt((wind.YGSE-wind.noseYGSE).^2 + (wind.ZGSE-wind.noseZGSE).^2) > 20 ...
%     & sqrt((ace.YGSE-ace.noseYGSE).^2 + (ace.ZGSE-ace.noseZGSE).^2) < 20; 

% When WIND is close to nose, we see the more accurate data (more SW
% effect)
% filter = sqrt((wind.YGSE-ace.YGSE).^2 + (wind.ZGSE-ace.ZGSE).^2) > 20 & ...
%         sqrt((wind.YGSE-wind.noseYGSE).^2 + (wind.ZGSE-wind.noseZGSE).^2) < 20;

% When GEOTAIL is close to nose, we see the more accurate data (more SW
% effect)
% filter = sqrt((geotail.YGSE-geotail.noseYGSE).^2 + (geotail.ZGSE-geotail.noseZGSE).^2) < 20;

% When WIND is close to nose, we see the more accurate data (more SW
% effect)
filter = sqrt((wind.YGSE-wind.noseYGSE).^2 + (wind.ZGSE-wind.noseZGSE).^2) > 20;

errorPercentageThreshold = 10000;
error = abs(ace.E_kl-wind.E_kl).*10^-3;
X1 = wind.E_kl.*10^-3;
errorPercentage = 100.*(error./X1);
errorPercentage(errorPercentage==Inf)=nan;

X2 = X1(errorPercentage<=errorPercentageThreshold & filter);
time2=datenum(wind.datetime(errorPercentage<=errorPercentageThreshold & filter));

% X2 = geotail.E_kl.*10^-3;
% time2=datenum(geotail.datetime);

% X2 = X1(error<=5);
% time2=datenum(wind.datetime(error<=5));


Y2 = omni.Fsml(time2);

[xindx2,E2]=discretize(X2,E);

for i = 1:max(xindx2)
        Y2gX2(i) = nanmean(Y2(xindx2==i));
        NSamples2(i) = sum(xindx2==i & ~isnan(Y2));
        SEM2(i) = nanstd(Y2(xindx2==i))./sqrt(NSamples2(i));
        ts2(i,:) = tinv([0.025 0.975], NSamples2(i)-1);
        CI2(i,:) = Y2gX2(i) + ts2(i,:)*SEM2(i);
        XBins2(i) = 0.5*(E2(i)+E2(i+1));
end

figure; 
CI1 = interp_nans(CI);
p1=plot(XBins, YgX,'k');
hold on;
plot_ci(XBins,CI1,'k',0.2);
ylim([-1800 0]);
xlim([0 30]);
hold on;


CI2 = interp_nans(CI2);
p2=plot(XBins2, Y2gX2,'r');
hold on;
plot_ci(XBins2,CI2,'r',0.2);
hold on;
plot(XBins2, smooth(CI2(:,1),10),'--m');

xlabel({'E_m [mV/m]','Geoeffective electric field'});
ylabel({'<SML|E_m> [nT]','Westward auroral electrojet current'});
legend([p1,p2],{'WIND','WIND >20RE'},'Location','southwest');
% legend(p1,{'1998-2019'},'Location','southwest');

%%
figure;
CI1 = interp_nans(CI);
p1=plot(XBins, YgX,'k');
hold on;
plot_ci(XBins,CI1,'k',0.2);
ylim([-1800 0]);
xlim([0 20]);

xlabel({'E_m [mV/m]','Geoeffective electric field'});
ylabel({'<SML|E_m> [nT]','Westward auroral electrojet current'});
legend(p1,{'1998-2019','Error in E_m < 5%'},'Location','southwest');
% hold on;
% 
% hold on;
% 
% CIy = interp_nans(CIy);
% p3=plot( XgY, YBins, 'b');
% hold on;
% plot_ci(interp_nans(XgY), CIy, 'b', 0.2);


%% See how the data is distributed between the different error percentages
figure; 
plot(2*ones(size(time2)),wind.datetime(errorPercentage<=errorPercentageThreshold),'.'); 
hold on; plot(ones(size(wind.datetime)),wind.datetime,'.'); 
xlim([0 3]);

%% Checking the ACE & WIND measurement in time series
time=wind.datetime;
indx = ~(isnan(ace.E_kl)) & ~(isnan(wind.E_kl)) & ~(isnan(geotail.E_kl)) &...
    time>datetime('16 Aug 1999 22:00','Format','dd MMM yyyy HH:mm') &...
    time<datetime('17 Aug 1999 12:00','Format','dd MMM yyyy HH:mm');
figure; 
subplot(2,1,1)
plot(time(indx),ace.E_kl(indx).*10^-3); 
hold on; 
plot(time(indx),wind.E_kl(indx).*10^-3);
hold on;
plot(time(indx),geotail.E_kl(indx).*10^-3);
legend('ACE','WIND','GEOTAIL');
ylim([0,5]);
ylabel('E_m');
subplot(2,1,2)
plot(time(indx), sqrt((ace.noseYGSE(indx)-ace.YGSE(indx)).^2 + (ace.noseZGSE(indx)-ace.ZGSE(indx)).^2)); hold on; 
plot(time(indx),sqrt((wind.noseYGSE(indx)-wind.YGSE(indx)).^2 + (wind.noseZGSE(indx)-wind.ZGSE(indx)).^2));
hold on;
plot(time(indx),sqrt((geotail.noseYGSE(indx)-geotail.YGSE(indx)).^2 + (geotail.noseZGSE(indx)-geotail.ZGSE(indx)).^2));
ylabel('Distance from Nose in Y-Z plane');
ylim([0,100]);


indx = ~(isnan(ace.E_kl)) & ~(isnan(wind.E_kl)) & ~(isnan(geotail.E_kl)) &...
    time>datetime('06 Aug 1999 22:00','Format','dd MMM yyyy HH:mm') &...
    time<datetime('07 Aug 1999 12:00','Format','dd MMM yyyy HH:mm');
figure; 
subplot(2,1,1)
plot(time(indx),ace.E_kl(indx).*10^-3); 
hold on; 
plot(time(indx),wind.E_kl(indx).*10^-3);
hold on;
plot(time(indx),geotail.E_kl(indx).*10^-3);
legend('ACE','WIND','GEOTAIL');
ylim([0,5]);
ylabel('E_m');
subplot(2,1,2)
plot(time(indx), sqrt((ace.noseYGSE(indx)-ace.YGSE(indx)).^2 + (ace.noseZGSE(indx)-ace.ZGSE(indx)).^2)); hold on; 
plot(time(indx),sqrt((wind.noseYGSE(indx)-wind.YGSE(indx)).^2 + (wind.noseZGSE(indx)-wind.ZGSE(indx)).^2));
hold on;
plot(time(indx),sqrt((geotail.noseYGSE(indx)-geotail.YGSE(indx)).^2 + (geotail.noseZGSE(indx)-geotail.ZGSE(indx)).^2));
ylabel('Distance from Nose in Y-Z plane');
ylim([0,100]);

%% Plot condfidence interval
function plot_ci(x,ci,color,alpha)
    hold on;
    X2 = [x, fliplr(x)];
    inBetween = [ci(:,1)', fliplr(ci(:,2)')];
    fill(X2,inBetween,color,'LineStyle','none','FaceAlpha',alpha);
end

% Extracting Omni data and developing gridded interpolants of the database
function omni = extract_omni_data(omniFile)
    
    %Time
    time = unixtime2matlab(h5read(omniFile,'/Time'));
    omni.time = time; 
    
        % Auroral electrojet indices
    SML = h5read(omniFile,'/Indices/SML');
    omni.Fsml = griddedInterpolant(time,SML);
    
    SMU = h5read(omniFile,'/Indices/SMU');
    omni.Fsmu = griddedInterpolant(time,SMU);
    
    AL = h5read(omniFile,'/Indices/AL');
    AL(AL==99999)=nan;
    omni.FAL = griddedInterpolant(time, AL);
    
    AU = h5read(omniFile,'/Indices/AU');
    AU(AU==99999)=nan;
    omni.FAU = griddedInterpolant(time, AU);
    
    symH = h5read(omniFile,'/Indices/SYM_H');
    omni.FsymH = griddedInterpolant(time, symH);
    % Solar wind
    
        % IMF Bz
    BzGSM = h5read(omniFile,'/BField/BzGSM');
    omni.FBz = griddedInterpolant(time,BzGSM);
    
        % IMF By
    ByGSM = h5read(omniFile,'/BField/ByGSM');
    omni.FBy = griddedInterpolant(time,ByGSM);
    
    
        % IMF B_T (Tangential to to GSM_x)
    BzGSM(BzGSM==0) = 0.0001;
    B_T = (ByGSM.^2 + BzGSM.^2).^0.5;
    
        % IMF Clock angle
    theta_kl = wrapTo2Pi(atan2(B_T,BzGSM));
        
        
        % Velocity
    velocity = h5read(omniFile,'/Velocity/V');
    omni.Fv = griddedInterpolant(time,velocity);
    
    
    % Solarwind - Magnetosphere Coupling
    l_0 = 7*(6371*10^3);
    
    E_kl = velocity.*B_T.*(sin(theta_kl/2)).^2; %Km nT/s
    % The “geoeffective” (or “merging”) electric field [Kan and Lee, 1979]
    omni.Fekl = griddedInterpolant(time,E_kl);    
    
end

function T = extract_set_of_data(loadFolder,satellite)
    fileStr = get_files_in_folder(loadFolder,[satellite,'_min_b*.txt']);
    for i=1:1:length(fileStr)
        T1 = extract_data([loadFolder,fileStr{i}]);
        if i==1
            T=T1;
        else
            T = [T;T1];
        end
    end
end


function T1 = extract_data(loadFile)
 
        format = '%4f %4f %3f %3f %4f %4f %4.1f %6f %6.2f %6.2f %6.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %6.2f %8.2f %8.2f %4f %8.1f %8.1f %8.1f %8.1f %7.2f %9.0f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %7f';
        tempData = load_ascii_files(loadFile, format, 0);
        T.datetime = datetime(tempData{1},1,tempData{2},tempData{3},tempData{4},zeros(size(tempData{4})));
        T.timeshift = tempData{8}; T.timeshift(T.timeshift==999999) = nan;
        T.BxGSM = tempData{13};  T.BxGSM(T.BxGSM > 9999) = nan;
        T.ByGSM = tempData{16}; T.ByGSM(T.ByGSM > 9999) = nan;
        T.BzGSM = tempData{17}; T.BzGSM(T.BzGSM > 9999) = nan;
        T.velocity = tempData{23}; T.velocity(T.velocity > 99999) = nan;
        T.B_T = (T.ByGSM.^2 + T.BzGSM.^2).^0.5;
        T.theta_kl = wrapTo2Pi(atan2(T.B_T,T.BzGSM));
        T.E_kl = T.velocity.*T.B_T.*(sin(T.theta_kl/2)).^2;
        T.XGSE = tempData{29}; T.XGSE(T.XGSE > 9999) = nan;
        T.YGSE = tempData{30}; T.YGSE(T.YGSE > 9999) = nan;
        T.ZGSE = tempData{31}; T.ZGSE(T.ZGSE > 9999) = nan;
        T.noseXGSE = tempData{32}; T.noseXGSE(T.noseXGSE > 9999) = nan;
        T.noseYGSE = tempData{33}; T.noseYGSE(T.noseYGSE > 9999) = nan;
        T.noseZGSE = tempData{34}; T.noseZGSE(T.noseZGSE > 9999) = nan;
        T.distance = sqrt((T.XGSE-T.noseXGSE).^2 + (T.YGSE-T.noseYGSE).^2 + (T.ZGSE-T.noseZGSE).^2);
        T1 = struct2table(T);
end

