%% Saturation of AE given E_M

%% Extracting data from ACE, WIND, and GEOTAIL
tic
dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
ace1 = extract_set_of_data(dataFolder,'ace');
wind1 = extract_set_of_data(dataFolder,'wind');
geotail1 = extract_set_of_data(dataFolder,'geotail');
omni = extract_omni_data(omniFile);
toc

%% ACE, WIND, GEOTAIL, time conjunctions
[~,~,wi]=intersect(geotail1.datetime,wind1.datetime);
[~,ai,wii]=intersect(ace1.datetime,wind1.datetime(wi));
[~,gi,~]=intersect(geotail1.datetime,wind1.datetime(wi(wii)));
wi = wi(wii);
%% ACE-WIND conjunction
[~,ai1,wi1]=intersect(ace1.datetime,wind1.datetime);
%% Geotail-Wind conjunction
[~,gi1,wi2] = intersect(geotail1.datetime,wind1.datetime);

%% ACE, WIND, GEOTAIL conjunction tables
ace = ace1(ai,:);
wind = wind1(wi,:);
geotail = geotail1(gi,:);

%% Filters
filter.ace = ((ace.noseYGSE-ace.YGSE).^2 + (ace.noseZGSE-ace.ZGSE).^2).^0.5 < 10;
filter.wind = ((wind.noseYGSE-wind.YGSE).^2 + (wind.noseZGSE-wind.ZGSE).^2).^0.5 < 10;
filter.geo = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 10;

filter.ace1 = ((ace1.noseYGSE-ace1.YGSE).^2 + (ace1.noseZGSE-ace1.ZGSE).^2).^0.5 < 10;
filter.wind1 = ((wind1.noseYGSE-wind1.YGSE).^2 + (wind1.noseZGSE-wind1.ZGSE).^2).^0.5 < 10;
filter.geo1 = ((geotail1.noseYGSE-geotail1.YGSE).^2 + (geotail1.noseZGSE-geotail1.ZGSE).^2).^0.5 < 10;


%% Plot all available data for Solar cycle 23

filter0 = ace1.datetime>datetime(1996,5,1) & ace1.datetime<datetime(2008,1,1) &...
    sqrt((ace1.YGSE-ace1.noseYGSE).^2 + (ace1.ZGSE-ace1.noseZGSE).^2) < 20;
curve1 = create_curve(ace1.E_kl(filter0).*10^-3,  omni.Fsml(datenum(ace1.datetime(filter0))), 30);

filter = wind1.datetime>datetime(1996,5,1) & wind1.datetime<datetime(2008,1,1) &...
    sqrt((wind1.YGSE-wind1.noseYGSE).^2 + (wind1.ZGSE-wind1.noseZGSE).^2) < 20;
time2=datenum(wind1.datetime(filter));
curve2 = create_curve(wind1.E_kl(filter).*10^-3, omni.Fsml(time2), curve1.E);

filter1 = geotail1.datetime>datetime(1996,5,1) & geotail1.datetime<datetime(2008,1,1) &...
    sqrt((geotail1.YGSE-geotail1.noseYGSE).^2 + (geotail1.ZGSE-geotail1.noseZGSE).^2) < 10;
time3=datenum(geotail1.datetime(filter1));
curve3 = create_curve(geotail1.E_kl(filter1).*10^-3, omni.Fsml(time3), curve1.E);

figure; 

subplot(2,2,1)

p1 = plot_curve(curve1,'k');
hold on;
p2 = plot_curve(curve2,'r');
hold on;
p3 = plot_curve(curve3,'b');
title('Solar Cycle: 23, May 1996 to Jan 2008'); 
xlabel({'E_m [mV/m]','Geoeffective electric field'});
ylabel({'<SML|E_m> [nT]','Westward auroral electrojet current'});
legend([p1,p2,p3],{'ACE','WIND','GEOTAIL'},'Location','southwest');
ylim([-1800 0]);
xlim([0 30]);

subplot(2,2,2)
histogram(ace1.datetime(filter0 & ~isnan(ace1.E_kl)),'FaceColor','k');
xlim([datetime(1996,5,1),datetime(2008,1,1)]); 
title('ACE');

subplot(2,2,3)

histogram(wind1.datetime(filter & ~isnan(wind1.E_kl)),'FaceColor','r');
xlim([datetime(1996,5,1),datetime(2008,1,1)]);
title('WIND');

subplot(2,2,4)

histogram(geotail1.datetime(filter1 & ~isnan(geotail1.E_kl)),'FaceColor','b');
xlim([datetime(1996,5,1),datetime(2008,1,1)]);
title('GEOTAIL');


%
% Plot orbits
figure;
subplot(2,2,1)
plot_orbit(ace1, filter0, 'k');
title('Solar Cycle 23');
legend('ACE');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');


subplot(2,2,2)
plot_orbit(wind1, filter, 'r');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');
legend('WIND');

subplot(2,2,3)
plot_orbit(geotail1, filter1, 'b');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');
legend('GEOTAIL');

%% Plot conjugate data for Solar cycle 23

filter0 = ace.datetime>datetime(1996,5,1) & ace.datetime<datetime(2008,1,1) &...
    sqrt((ace.YGSE-ace.noseYGSE).^2 + (ace.ZGSE-ace.noseZGSE).^2) < 20;
curve1 = create_curve(ace.E_kl(filter0).*10^-3,  omni.Fsml(datenum(ace.datetime(filter0))), 30);

filter = wind.datetime>datetime(1996,5,1) & wind.datetime<datetime(2008,1,1) &...
    sqrt((wind.YGSE-wind.noseYGSE).^2 + (wind.ZGSE-wind.noseZGSE).^2) < 20;
time2=datenum(wind.datetime(filter));
curve2 = create_curve(wind.E_kl(filter).*10^-3, omni.Fsml(time2), curve1.E);

filter1 = geotail.datetime>datetime(1996,5,1) & geotail.datetime<datetime(2008,1,1) &...
    sqrt((geotail.YGSE-geotail.noseYGSE).^2 + (geotail.ZGSE-geotail.noseZGSE).^2) < 20;
time3=datenum(geotail.datetime(filter1));
curve3 = create_curve(geotail.E_kl(filter1).*10^-3, omni.Fsml(time3), curve1.E);

figure; 

subplot(2,2,1)

p1 = plot_curve(curve1,'k');
hold on;
p2 = plot_curve(curve2,'r');
hold on;
p3 = plot_curve(curve3,'b');
title('Solar Cycle: 23, May 1996 to Jan 2008'); 
xlabel({'E_m [mV/m]','Geoeffective electric field'});
ylabel({'<SML|E_m> [nT]','Westward auroral electrojet current'});
legend([p1,p2,p3],{'ACE','WIND','GEOTAIL'},'Location','southwest');
ylim([-1800 0]);
xlim([0 30]);

subplot(2,2,2)
histogram(ace.datetime(filter0 & ~isnan(ace.E_kl)),'FaceColor','k');
xlim([datetime(1996,5,1),datetime(2008,1,1)]); 
title('ACE');

subplot(2,2,3)

histogram(wind.datetime(filter & ~isnan(wind.E_kl)),'FaceColor','r');
xlim([datetime(1996,5,1),datetime(2008,1,1)]);
title('WIND');

subplot(2,2,4)

histogram(geotail.datetime(filter1 & ~isnan(geotail.E_kl)),'FaceColor','b');
xlim([datetime(1996,5,1),datetime(2008,1,1)]);
title('GEOTAIL');


%%
% Plot orbits
figure;
subplot(2,2,1)
plot_orbit(ace, filter0, 'k');
title('Solar Cycle 23');
legend('ACE');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');


subplot(2,2,2)
plot_orbit(wind, filter, 'r');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');
legend('WIND');

subplot(2,2,3)
plot_orbit(geotail, filter1, 'b');
xlim([-100,100]);
ylim([-50,50]);
xlabel('YGSE-YGSE_{nose}');
ylabel('ZGSE-ZGSE_{nose}');
legend('GEOTAIL');



%% 
% Reponse of SML to E_m from OMNI database
% filter = sqrt((wind1.YGSE-wind1.noseYGSE).^2 + (wind1.ZGSE-wind1.noseZGSE).^2)<10; 
% Y = omni.Fsml(datenum(wind1.datetime));
% X = wind1.E_kl.*10^-3;
filter0 = ace1.datetime>datetime(2008,1,1) & ace1.datetime<datetime(2017,1,1);
Y = omni.Fsml(datenum(ace1.datetime(filter0)));
X = ace1.E_kl(filter0).*10^-3;

% Y = omni.Fsml(datenum(geotail1.datetime));
% X = geotail1.E_kl.*10^-3;
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

% filter = sqrt((geotail1.YGSE-geotail1.noseYGSE).^2 + (geotail1.ZGSE-geotail1.noseZGSE).^2) < 5;

% When WIND is close to nose, we see the more accurate data (more SW
% effect)
% filter = sqrt((wind1.YGSE-wind1.noseYGSE).^2 + (wind1.ZGSE-wind1.noseZGSE).^2) < 30 & wind1.datetime>datetime(1997,1,1);

filter = wind1.datetime>datetime(2008,1,1) & wind1.datetime<datetime(2017,1,1);

% errorPercentageThreshold = 10000;
% error = abs(ace1.E_kl-wind1.E_kl).*10^-3;
X1 = wind1.E_kl.*10^-3;
% X1 = geotail1.E_kl.*10^-3;
% errorPercentage = 100.*(error./X1);
% errorPercentage(errorPercentage==Inf)=nan;

X2 = X1(filter);
time2=datenum(wind1.datetime(filter));

% X2 = geotail.E_kl.*10^-3;
% time2=datenum(geotail1.datetime(filter));

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
legend([p1,p2],{'GEOTAIL','GEOTAIL <5RE'},'Location','southwest');

disp('Done');
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

function curve = create_curve(X, Y, Ei)
% Plot Y|X
if nargin<3
    Ei = 100;
end

Y(Y==999999)=nan;
[xindx, E] = discretize(X,Ei);

for i = 1:max(xindx)
    curve.YgX(i) = nanmean(Y(xindx==i));
    curve.NSamples(i) = sum(xindx==i & ~isnan(Y));
    curve.SEM(i) = nanstd(Y(xindx==i))./sqrt(curve.NSamples(i));
    curve.ts(i,:) = tinv([0.025 0.975],curve.NSamples(i)-1);
    curve.CI(i,:) = curve.YgX(i) + curve.ts(i,:)*curve.SEM(i);
    curve.XBins(i) = 0.5*(E(i)+E(i+1));
end

curve.E = E; 

end

function p = plot_curve(curve, color)

CI1 = interp_nans(curve.CI);
p=plot(curve.XBins, curve.YgX,color);
hold on;
plot_ci(curve.XBins,CI1,color,0.2);
end

function p = plot_orbit(satellite, filter, color)

p = scatter(satellite.noseYGSE(filter)-satellite.YGSE(filter),...
    satellite.noseZGSE(filter)-satellite.ZGSE(filter),5,...
    'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.01);
end

function p = plot_orbit_2(satellite, color)

p = scatter(satellite.noseXGSE-satellite.XGSE,satellite.noseZGSE-satellite.ZGSE,5,...
    'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.01);
end
