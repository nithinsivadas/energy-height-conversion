%% Modeling influence of probability distribution on error


dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
geotail1 = extract_set_of_data(dataFolder,'geotail');
wind1 = extract_set_of_data(dataFolder,'wind');
%%
pdwind = fitdist(wind1.E_kl(~isnan(wind1.E_kl) & wind1.E_kl>0)*10^-3,'kernel');
pdgeo = fitdist(geotail1.E_kl(~isnan(geotail1.E_kl) & geotail1.E_kl>0)*10^-3,'kernel');
%%
Nsamples = 105120000;
Nensemble = 1;
pd = makedist('Lognormal',-0.286163,0.789428);
xTrue = random(pd,Nsamples,Nensemble); 
R = 2;
yTrue = 2.*xTrue; 


% xMeasured = xTrue + xTrue.*random('Normal',+0.08,0.25,Nsamples,Nensemble);
xMeasured = xTrue.*(random('Lognormal',0,0.25,Nsamples,Nensemble));
pdTrue = fitdist(xTrue(xTrue>0),'Lognormal');
pdMeasured = fitdist(xMeasured(xMeasured>0),'Lognormal');


%% Histogram
xBins = 0:0.5:30;
figure; 
histogram(xTrue,xBins,'Normalization','pdf');
hold on;
histogram(xMeasured(xMeasured>0),xBins,'Normalization','pdf');

hold on;
plot(xBins, pdf(pdTrue,xBins),'b','LineWidth',1);
hold on;
plot(xBins, pdf(pdMeasured,xBins),'r','LineWidth',1);

legend(['True: exp(-',num2str(pdTrue.mean,3),') X'],['Measured: exp(-',num2str(pdMeasured.mean,3),') X']);
xlabel('X');
ylabel('PDF');
%% Output calculation
yMeasured = yTrue + yTrue.*random('Normal',0,1,Nsamples,1);

%% Underestimation of conditional expectation 
curve1 = create_curve(xMeasured,yMeasured,xBins);
curve2 = create_curve(xTrue,yMeasured,xBins);
curve3 = create_curve(xTrue,yTrue,xBins);
curve11 = create_curve(xMeasured,yTrue,xBins);

figure;

pp1=plot_curve(curve1,'k');
hold on;

pp2=plot_curve(curve2,'r');
hold on;

pp3=plot_curve(curve3,'b');
hold on;

pp4=plot_curve(curve11,'g');
ylim([0,100]);
xlim([0,50]);
xlabel('X'); ylabel('Y');
legend([pp1,pp2,pp3,pp4],{'<Y_{Measured}|X_{Measured}>','<Y_{Measured}|X_{True}>','<Y_{True}|X_{True}>','<Y_{True}|X_{Measured}>'});

%% Underestimating the true input, due to spill over from exponential distribution
curve4 = create_curve(xMeasured,xTrue,xBins);
curve5 = create_curve(xMeasured,xMeasured,xBins);

%%
figure; 

p1=plot_curve(curve4,'k');
hold on;

p2=plot_curve(curve5,'b');
ylim([0,30]);
xlim([0,30]);
legend([p1,p2],{'<X_{True}|X_{Measured}>','<X_{Measured}|X_{Measured}>'},'Location','northwest');
xlabel('X_{Measured}');

%%
figure; 
histogram(xTrue(xMeasured<10.5 & xMeasured>9.5),0:0.5:30,'Normalization','pdf');
hold on;
histogram(xMeasured(xMeasured<10.5 & xMeasured>9.5),0:0.5:30,'Normalization','pdf');
legend('True','Measured');
xlabel('X');
ylabel('pdf');

%%
dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
% ace1 = extract_set_of_data(dataFolder,'ace');
wind1 = extract_set_of_data(dataFolder,'wind');
geotail1 = extract_set_of_data(dataFolder,'geotail');
%% Probability distribution of measurements

figure; 
histogram(geotail1.E_kl(~isnan(geotail1.E_kl))*10^-3,xBins,'Normalization','pdf');
hold on;
histogram(wind1.E_kl(~isnan(wind1.E_kl))*10^-3,xBins,'Normalization','pdf');

pdgeo = fitdist(geotail1.E_kl(~isnan(geotail1.E_kl) & geotail1.E_kl>0)*10^-3,'Gamma');
pdwind = fitdist(wind1.E_kl(~isnan(wind1.E_kl) & wind1.E_kl>0)*10^-3,'Gamma');
hold on;
plot(xBins,pdf(pdgeo,xBins),'b','LineWidth',1);
hold on;
plot(xBins,pdf(pdwind,xBins),'r','LineWidth',1);
ylim([0,1]);
legend('Geotail (True)','Wind (Measured)',...
    ['Geotail: A : ',num2str(pdgeo.a,3),', B : ',num2str(pdgeo.b,3)],...
    ['Wind: A : ',num2str(pdwind.a,3),', B : ',num2str(pdwind.b,3)]);

xlabel('E_{kl} [mV/m]');
ylabel('pdf');

%% Finding the probability distribution of solar wind electric field
%% Geotail1
time = geotail1.datetime; 
E_kl_1min = geotail1.E_kl.*10^-3;
distance_1min = ((geotail1.noseYGSE-geotail1.YGSE).^2 + (geotail1.noseZGSE-geotail1.ZGSE).^2).^0.5;

[time_1hr,E_kl_1hr] = hourly_bin(time,E_kl_1min);
[~,distance_1hr] = hourly_bin(time,distance_1min);


figure; 
plot(time_1hr(distance_1hr<7),E_kl_1hr(distance_1hr<7));
hold on;
plot(time(distance_1min<7),E_kl_1min(distance_1min<7));

%% Plot histogram
% XBins = logspace(-3,+3,100);
XBins = 0:0.5:30;
y1 = E_kl_1min(~isnan(E_kl_1min) & distance_1min<10 & E_kl_1min>0 );
y2 = E_kl_1hr(~isnan(E_kl_1hr) & distance_1hr<10 & E_kl_1hr>0);

pdf_1min = fitdist(y1,'Lognormal');
pdf_1hr = fitdist(y2,'Lognormal'); 
%Note: A lognormal fits for 1.5 hr binned geotail data, distances<10RE from nose, 
% mu = -0.286163 ; sigma = 0.789428  [0.776011, 0.803321]

% [h_1min,p_1min,stat1] = chi2gof(y1,'CDF',pdf_1min);
% [h_1hr,p_1hr,stat2] = chi2gof(y2,'CDF',pdf_1hr);

% [h_1min,p_1min,k_1min] = kstest2(y1,random(pdf_1min,10000,1))
[h_1min,p_1min,k_1min] = kstest(y1,'CDF',pdf_1min)
% [h_1hr,p_1hr,k_1hr] = kstest2(y2,random(pdf_1hr,100000,1).*(1 + random(gmdistribution(0,0.25),100000)))
[h_1hr,p_1hr,k_1hr] = kstest(y2,'CDF',pdf_1hr)

figure;
histogram(y1,XBins,'Normalization','pdf');
hold on;
histogram(y2,XBins,'Normalization','pdf');
hold on;
plot(XBins,pdf(pdf_1min,XBins),'b');
hold on;
plot(XBins,pdf(pdf_1hr,XBins),'r');

set(gca,'XScale','linear');
xlabel('E_{kl} [mV/m]');
ylabel('pdf');
legend('Geotail, 1min Bins, <11 RE from nose',...
    'Geotail, 1.5hr Bins, <11 RE from nose',...
    ['1-min exp fit, Fail chi-2, \mu = ',num2str(pdf_1min.mean,3)],...
    ['1.5-hr exp fit, Pass chi-2, \mu = ',num2str(pdf_1hr.mean,3)]);




%% Wind
time = wind1.datetime; 
E_kl_1min = wind1.E_kl.*10^-3;
distance_1min = ((wind1.noseYGSE-wind1.YGSE).^2 + (wind1.noseZGSE-wind1.ZGSE).^2).^0.5;

[time_1hr,E_kl_1hr] = hourly_bin(time,E_kl_1min);
[~,distance_1hr] = hourly_bin(time,distance_1min);


figure; 
plot(time_1hr(distance_1hr<7),E_kl_1hr(distance_1hr<7));
hold on;
plot(time(distance_1min<7),E_kl_1min(distance_1min<7));

%% Plot histogram
XBins = 0:0.5:50;
y1 = E_kl_1min(~isnan(E_kl_1min) & distance_1min<1000 & E_kl_1min>0 );
y2 = E_kl_1hr(~isnan(E_kl_1hr) & distance_1hr<1000 & E_kl_1hr>0);

pdf_1min = fitdist(y1,'Lognormal');
pdf_1hr = fitdist(y2,'Lognormal');

[h_1min,p_1min,k_1min] = kstest(y1,'CDF',pdf_1min)
[h_1hr,p_1hr,k_1hr] = kstest(y2,'CDF',pdf_1hr)

figure;
histogram(y1,XBins,'Normalization','pdf');
hold on;
histogram(y2,XBins,'Normalization','pdf');
hold on;
plot(XBins,pdf(pdf_1min,XBins),'b');
hold on;
plot(XBins,pdf(pdf_1hr,XBins),'r');

%%
pp = fitdist(E_kl_1hr(~isnan(E_kl_1hr))*10^-3,'Weibull');
h = chi2gof(E_kl_1hr(~isnan(E_kl_1hr))*10^-3,'CDF',pp)
%%
hold on; 
plot(1:1:30,pdf(pp,1:1:30));

function curve = create_curve(X, Y, Ei)

if nargin<3
    Ei = 100;
end
Y = Y(:);
Y(Y==999999)=nan;
[xindx, E] = discretize(X(:),Ei);

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

function plot_ci(x,ci,color,alpha)
    hold on;
    X2 = [x, fliplr(x)];
    inBetween = [ci(:,1)', fliplr(ci(:,2)')];
    fill(X2,inBetween,color,'LineStyle','none','FaceAlpha',alpha);
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

function [X_1hr,Y_1hr] = hourly_bin(X,Y)
    [~,~,idx] = histcounts(X,X(1):hours(1.5):X(end-1));
    Y_1hr = accumarray(idx(idx~=0),Y(idx~=0),[],@(x)mean(x,'omitnan'));
    X_1hr = accumarray(idx(idx~=0),datenum(X(idx~=0)),[],@(x)mean(x,'omitnan'));
    X_1hr = datetime(X_1hr, 'ConvertFrom','datenum');
end