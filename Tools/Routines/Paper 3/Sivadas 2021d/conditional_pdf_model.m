%% Paper 6: P(X|Xm)
%% Initializing some data
dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
geotail1 = extract_set_of_data(dataFolder,'geotail');
wind1 = extract_set_of_data(dataFolder,'wind');
ace1 = extract_set_of_data(dataFolder,'ace');

%%
Nsamples = 105120000;
Nensemble = 1;

% Creating the random variable X
pdfx = makedist('Lognormal',-0.286163,0.789428);
% pdfx = makedist('Uniform',0,30);
% pdfx = makedist('Exponential',1.04);
X = random(pdfx,Nsamples,Nensemble); % Samples/instances of the random variable
sig=0.15;

% Reconstructing the pdf of random variable X from its samples
% pdfx1 = fitdist(X(X>0),'Lognormal');

% Inducing measurement error Type3

% Two separate, uncorrelated measurements of X
% Xm_3a = X + random('Normal',+0.07,0.12.*X,Nsamples,Nensemble); 
% Xm_3b = X + random('Normal',+0.07,0.12.*X,Nsamples,Nensemble);
% Xm_3a = X + random('Normal',+0.07,0.24.*X,Nsamples,Nensemble); 
% Xm_3b = X + random('Normal',+0.07,0.24.*X,Nsamples,Nensemble);

Xm_3a = X + X.*random(pdf1,Nsamples,Nensemble); 
Xm_3b = X + X.*random(pdf1,Nsamples,Nensemble); 

% Xm_3a = X + random('tlocationscale',0,X.*sig./2,2,Nsamples,Nensemble);
% Xm_3b = X + random('tlocationscale',0,X.*sig./2,2,Nsamples,Nensemble);

% Xm_3a = X.*exp(random('tlocationscale',log(1.05./(1+sig.^2).^0.5),sig./2,1.31,Nsamples,Nensemble));
% Xm_3b = X.*exp(random('tlocationscale',log(1.05./(1+sig.^2).^0.5),sig./2,1.31,Nsamples,Nensemble));

% Xm_3a = X.*random('tlocationscale',1,sig./2,1.31,Nsamples,Nensemble);
% Xm_3b = X.*random('tlocationscale',1,sig./2,1.31,Nsamples,Nensemble);

% Xm_3a = X.*random('Lognormal',0,0.5,Nsamples,Nensemble);
% Xm_3b = X.*random('Lognormal',0,0.5,Nsamples,Nensemble);
% pdfm3 = fitdist(Xm3(Xm3>0),'Lognormal'); % Reconstructing its pdf, not the right pdf

% Plotting P(e|Xm)

% Estimating the error
e1 = Xm_3a - Xm_3b; %Type 2 - difference between two erroneous measurements
e = Xm_3a - X; %Type 1 - difference between erroneous and true-value


[fe1gXm_3a,XX,YY] = conditional_pdf(e1,Xm_3a,-5:0.1:5,0:1:30);
[fegXm_3a,XX1,YY1] = conditional_pdf(e,Xm_3a,-5:0.1:5,0:1:30);

[fe1gX,XXT,YYT] = conditional_pdf(e1,X,-5:0.1:5,0:1:30);
[fegX,XXT1,YYT1] = conditional_pdf(e,X,-5:0.1:5,0:1:30);
%
% [ee3,XXm3,cpdf_egxm3] = conditional_pdf(e,Xm3,-5:0.1:5,0:1:20);
h=figure; 
resize_figure(h,200,400);

subplot(2,4,1)
p=pcolor(XX,YY,fe1gXm_3a);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X*');
xlabel('P(X*_1-X*_2|X*_1)');

subplot(2,4,2)
p=pcolor(XXT,YYT,fe1gX);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X');
xlabel('P(X*_1-X*_2|X)');

subplot(2,4,3)
p=pcolor(XX1,YY1,fegXm_3a);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X*');
xlabel('P(X*_1-X|X*_1)');

subplot(2,4,4)
p=pcolor(XXT1,YYT1,fegX);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X');
xlabel('P(X*_1-X|X)');

% Plotting distribution 
subplot(2,4,5)
histogram(X,logspace(-3,3,100),'Normalization','pdf'); 
set(gca,'XScale','log');
hold on;
histogram(Xm_3a,logspace(-3,3,100),'Normalization','pdf'); 
legend('X','X*');

% Calculating conditional expectations
xBins = 0:0.5:30;
curve31 = create_curve(Xm_3a,X,xBins);
curve32 = create_curve(X,X,xBins);

% Plotting conditional expectation 
subplot(2,4,6)
p1=plot_curve(curve32,'k');
hold on;
p2=plot_curve(curve31,'b');
ylim([0,30]);
xlim([0,30]);
legend([p1,p2],{'<X|X>','<X|X*>'},'Location','northwest');
xlabel('X*');
ylabel('Conditional expectation of X');

% [fXgXm_3a,XXm1,YYm1] = conditional_pdf(X,Xm_3a,0:1:30,0:1:30);
% 
% subplot(2,4,7)
% p=pcolor(XXm1,YYm1,fXgXm_3a');
% set(p, 'EdgeColor','none');
% colorbar;
% caxis([0,0.5]);
% xlabel('X*');
% ylabel('P(X|X*)');
% hold on;
% plot(0:0.1:30,0:0.1:30,'w','LineWidth',2);
% colormap(inferno);

R = 2; 
Y = R.*X + random('Normal',0,1,Nsamples,Nensemble);
Ym_3a = R.*Xm_3a + random('Normal',0,1,Nsamples,Nensemble);

% Calculating conditional expectations for Y|X
curve31y = create_curve(Xm_3a,Y,xBins);
curve32y = create_curve(X,Y,xBins);

% Plotting conditional expectation 
subplot(2,4,7)
p3=plot_curve(curve32y,'k');
hold on;
p4=plot_curve(curve31y,'b');
ylim([0,60]);
xlim([0,30]);
legend([p3,p4],{'<Y|X>','<Y*|X*>'},'Location','northwest');
xlabel('X*');
ylabel('Conditional expectation of Y|X');

[fYgX,YYY,XXX] = conditional_pdf(Y,X,0:1:20,0:1:20);
[fYgXm_3a,YYYm,XXXm] = conditional_pdf(Y,Xm_3a,0:1:30,0:1:30);

subplot(2,4,8)
p=pcolor(YYYm,XXXm,fYgXm_3a');
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.5]);
xlabel('X*');
ylabel('P(Y|X*)');
hold on;
plot(0:0.1:30,R.*(0:0.1:30),'w','LineWidth',2);
colormap(inferno);


%% Checkng the value of error from WIND, ACE, GEOTAIL 

%% ACE, WIND, GEOTAIL, time conjunctions
[~,~,wi]=intersect(geotail1.datetime,wind1.datetime);
[~,ai,wii]=intersect(ace1.datetime,wind1.datetime(wi));
[~,gi,~]=intersect(geotail1.datetime,wind1.datetime(wi(wii)));
wi = wi(wii);

% ACE-WIND conjunction
[~,ai1,wi1]=intersect(ace1.datetime,wind1.datetime);
% Geotail-Wind conjunction
[~,gi1,wi2] = intersect(geotail1.datetime,wind1.datetime);

% ACE, WIND, GEOTAIL conjunction tables
ace = ace1(ai,:);
wind = wind1(wi,:);
geotail = geotail1(gi,:);
%%
% Filters
filter.ace = ((ace.noseYGSE-ace.YGSE).^2 + (ace.noseZGSE-ace.ZGSE).^2).^0.5 < 10;
filter.wind = ((wind.noseYGSE-wind.YGSE).^2 + (wind.noseZGSE-wind.ZGSE).^2).^0.5 < 10;
filter.geo = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 7;

filter.ace1 = ((ace1.noseYGSE-ace1.YGSE).^2 + (ace1.noseZGSE-ace1.ZGSE).^2).^0.5 < 10;
filter.wind1 = ((wind1.noseYGSE-wind1.YGSE).^2 + (wind1.noseZGSE-wind1.ZGSE).^2).^0.5 < 10;
filter.geo1 = ((geotail1.noseYGSE-geotail1.YGSE).^2 + (geotail1.noseZGSE-geotail1.ZGSE).^2).^0.5 < 10;

%
filter.geo = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 20;
XmS_a = wind.E_kl(filter.geo)*10^-3;
XmS_b = ace.E_kl(filter.geo)*10^-3;
XS = geotail.E_kl(filter.geo)*10^-3;
% XmS_a = XS + random('Normal',0,0.25.*XS,length(XS),1); 
% XmS_b = XS + random('Normal',0,0.25.*XS,length(XS),1); 
e1S = XmS_a-XmS_b;
eS = XmS_a - XS; 

[fe1gXmS_a,XXS,YYS] = conditional_pdf(e1S,XmS_a,-5:0.1:5,0:1:20);

[fegXmS_a,XX1S,YY1S] = conditional_pdf(eS,XmS_a,-5:0.1:5,0:1:20);


[fe1gXS,XXTS,YYTS] = conditional_pdf(e1S,XS,-5:0.1:5,0:1:20);
[fegXS,XXT1S,YYT1S] = conditional_pdf(eS,XS,-5:0.1:5,0:1:20);

figure;
subplot(2,4,1)
p=pcolor(XXS,YYS,fe1gXmS_a);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X*_{wind}');
xlabel('P(X*_{wind}-X*_{ace}|X*_{wind})');
colormap(inferno);


subplot(2,4,2)
p=pcolor(XX1S,YY1S,fegXmS_a);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X*_{wind}');
xlabel('P(X*_{wind}-X_{geotail}|X*_{wind})');
colormap(inferno);

subplot(2,4,3)
p=pcolor(XXTS,YYTS,fe1gXS);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P(X*_{wind}-X*_{ace}|X_{geotail})');
colormap(inferno);


subplot(2,4,4)
p=pcolor(XXT1S,YYT1S,fegXS);
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P(X*_{wind}-X_{geotail}|X_{geotail})');
colormap(inferno);

% Plotting distribution 
subplot(2,4,5)
histogram(XS,logspace(-3,3,100),'Normalization','pdf'); 
set(gca,'XScale','log');
hold on;
histogram(XmS_a,logspace(-3,3,100),'Normalization','pdf'); 
legend('X_{geotail}','X*_{wind}');

% Calculating conditional expectations
curve31S = create_curve(XmS_a,XS,xBins);
curve32S = create_curve(XS,XS,xBins);

% Plotting conditional expectation 
subplot(2,4,6)
p1=plot_curve(curve32S,'k');
hold on;
p2=plot_curve(curve31S,'b');
ylim([0,30]);
xlim([0,30]);
legend([p1,p2],{'<X_{geotail}|X_{geotail}>','<X_{geotail}|X*_{wind}>'},'Location','northwest');
xlabel('X');
ylabel('Conditional expectation of X');

[fXgXmS_a,XXm1S,YYm1S] = conditional_pdf(XS,XmS_a,0:1:20,0:1:20);

subplot(2,4,7)
p=pcolor(XXm1S,YYm1S,fXgXmS_a');
set(p, 'EdgeColor','none');
colorbar;
caxis([0,1]);
xlabel('X*_{wind}');
ylabel('P(X_{geotail}|X*_{wind})');
hold on;
plot(0:0.1:20,0:0.1:20,'w','LineWidth',2);
colormap(inferno);

R = 2; 
YS = R.*XS + random('Normal',0,1,length(XS),1);
YmS_a = R.*XmS_a + random('Normal',0,1,length(XS),1);

[fYSgXS,YYYS,XXXS] = conditional_pdf(YS,XS,0:1:20,0:1:20);
[fYSgXmS_a,YYYmS,XXXmS] = conditional_pdf(YS,XmS_a,0:1:20,0:1:20);

subplot(2,4,8)
p=pcolor(YYYmS,XXXmS,fYSgXmS_a');
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.5]);
xlabel('X*_{wind}');
ylabel('P(Y_{SML}|X*_{wind})');
hold on;
plot(0:0.1:20,R.*(0:0.1:20),'w','LineWidth',2);
colormap(inferno);

%%
figure;

p=pcolor(YYYmS,XXXmS,fYSgXmS_a');
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.1]);
xlabel('X*_{wind}');
ylabel('P(Y_{SML}|X*_{wind})');
hold on;
plot(0:0.1:20,R.*(0:0.1:20),'w','LineWidth',2);
colormap(inferno);


%% Error estimate from data
i=200; f1 = fit(XXT1S(i,:).',fegXS(i,:).','gauss1'); figure; plot(XXT1S(i,:),(fegXS(i,:))); hold on; plot(XXT1S(i,:),f1(XXT1S(i,:))); title(['X = ',num2str(YYT1S(i,1)),'; \sigma = ' num2str(f1.c1),'; Ratio =',num2str(f1.c1./YYT1S(i,1))]);

figure; plot((std(fegXmS_a')./YYT1S(:,1)')); ylim([0,1]);

%% Calculating sigma_TM
% Error in propagation of L1 readings to magnetopause nose

time = geotail1.datetime(100:end-100);
X4 = geotail1.E_kl(100:end-100).*10^-3;
sigmaMT = 1;
timeMT = time + minutes(random('Normal',0,sigmaMT,length(time),1));

sigmaTI = 20;
timeTI = time + minutes(random('Normal',5,sigmaTI,length(time),1));

X4MT = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeMT));
eMT = X4-X4MT;

X4TI = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeTI));
eTI = X4-X4TI; 

timeT = time + minutes(random('Normal',0,sigmaMT,length(time),1)) +...
    minutes(random('Normal',0,sigmaTI,length(time),1));
X4T = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeT));
eT = X4-X4T;

[feMTgXS,XMTS,YMTS] = conditional_pdf(eMT,X4,-20:0.5:20,0:1:20);
[feTIgXS,XTIS,YTIS] = conditional_pdf(eTI,X4,-20:0.5:20,0:1:20);
[feTgXS,XTS,YTS] = conditional_pdf(eT,X4,-20:0.5:20,0:1:20);

figure;
p=pcolor(XMTS,YMTS,feMTgXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{TM}|X_{geotail})');
colormap(inferno);
title('\sigma_{TM} = 2 min');
xlim([-5,5]);

figure;
p=pcolor(XTIS,YTIS,feTIgXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{TI}|X_{geotail})');
colormap(inferno);
title('\sigma_{TI} = 20 min');
xlim([-5,5]);

figure;
p=pcolor(XTS,YTS,feTgXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{T}|X_{geotail})');
title('\sigma_{T} = 2 + 20 min');
colormap(inferno);
xlim([-5,5]);

%% Plot std.dev of error varying with mean of the variable
curve41 = create_curve(X4, eMT, 0:0.5:30);
curve42 = create_curve(X4, eTI, 0:0.5:30);

figure; 
plot(curve41.XBins,100.*curve41.stdYgX./curve41.XBins); 
hold on;
plot(curve42.XBins,100.*curve42.stdYgX./curve42.XBins);
legend('\sigma_{TM}','\sigma_{TI}');
xlabel('X_{geotail}');
ylabel('\sigma /<X_{geotail}> [%]');

%% Calculating the std.dev of error (for MI, TI, and spatial error)
time = geotail.datetime(20:end-20);
X4 = geotail.E_kl(20:end-20).*10^-3;
X4a = ace.E_kl(20:end-20).*10^-3;
X4b = wind.E_kl(20:end-20).*10^-3;
sigmaMT = 1;
timeMT = time + minutes(random('Weibull',sigmaMT,2,length(time),1));

sigmaTI = 20;
timeTI = time + minutes(random('Weibull',sigmaTI,2,length(time),1));

X4MT = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeMT));
eMT = X4-X4MT;

X4TI = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeTI));
eTI = X4-X4TI; 

timeT = time + minutes(random('Weibull',sigmaMT,2,length(time),1)) +...
    minutes(random('Weibull',sigmaTI,2,length(time),1));
X4T = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeT));
eT = X4-X4T;
eS = X4 - X4b;
e = eT + eS;

curve41 = create_curve(X4, eMT, 0:0.5:30);
curve42 = create_curve(X4, eTI, 0:0.5:30);
curve43 = create_curve(X4, eS, 0:0.5:30);
curve44 = create_curve(X4, e, 0:0.5:30);

figure; 
plot(curve41.XBins,100.*curve41.stdYgX./curve41.XBins); 
hold on;
plot(curve42.XBins,100.*curve42.stdYgX./curve42.XBins);
hold on;
plot(curve43.XBins,100.*curve43.stdYgX./curve43.XBins);
hold on;
plot(curve44.XBins,100.*curve44.stdYgX./curve44.XBins,'k');
legend('Time to Magnetopause nose \sigma_{TM}','Time to Ionosphere \sigma_{TI}',...
    'Spatial variability \sigma_{Wind-Geotail}','Total Error \sigma_e');
ylim([0,100]);
xlabel({'Mean of approx. True Value','X_{geotail}'});
ylabel({'Error Percentage','\sigma /<X_{geotail}> [%]'});

%%
[fegXS,XeTS,YeTS] = conditional_pdf(e,X4,-20:0.5:20,0:1:20);
figure;
p=pcolor(XeTS,YeTS,fegXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.2]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{Total}|X_{geotail})');
title('\sigma_{T} = 2 + 20 min & \sigma_S = spatial error (WIND-GEOTAIL]');
colormap(inferno);
xlim([-5,5]);

[fegXS,XeTS,YeTS] = conditional_pdf(eT,X4,-20:0.5:20,0:1:20);
figure;
p=pcolor(XeTS,YeTS,fegXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.2]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{TI+TM}|X_{geotail})');
title('\sigma_{T} = 2 + 20 min - weibull distribution');
colormap(inferno);
xlim([-5,5]);

[fegXS,XeTS,YeTS] = conditional_pdf(eS,X4,-20:0.5:20,0:1:20);
figure;
p=pcolor(XeTS,YeTS,fegXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,0.2]);
ylabel('X_{geotail}');
xlabel('P(\epsilon_{X*-X}|X_{geotail})');
title('\sigma_S = spatial error [WIND-GEOTAIL]');
colormap(inferno);
xlim([-5,5]);

%% Compare the error distribution between measured and toy-model

% Estimating spacecraft error
time = geotail.datetime(20:end-20);
X4 = geotail.E_kl(20:end-20).*10^-3;
X4a = ace.E_kl(20:end-20).*10^-3;
X4b = wind.E_kl(20:end-20).*10^-3;
eS = X4 - X4b;

% Estimating toy-model error
Nsamples = length(time);
Nensemble = 1;

% Creating the random variable X
pdfx = makedist('Lognormal',-0.286163,0.789428);
X = random(pdfx,Nsamples,Nensemble); % Samples/instances of the random variable
sig = 0.15;
Xm_3a = X + X.*random('Normal',0,sig,Nsamples,Nensemble);
% Xm_3b = X.*random('Lognormal',log(1./(1+sig.^2).^0.5),(log(1+sig.^2)).^0.5,Nsamples,Nensemble);
Xm_3b = X + X.*random('tlocationscale',+0.01,sig./2,2,Nsamples,Nensemble);
% Xm_3b = X.*exp(random('tlocationscale',log(1.05./(1+sig.^2).^0.5),sig./2,0.5,Nsamples,Nensemble));
% Estimating the error
e = X-Xm_3b; %Type 1 - difference between erroneous and true-value
e1 = X-Xm_3a;

%%
sigmaMT = 8;
timeMT = time + minutes(random('Weibull',sigmaMT,2,length(time),1));

figure;
histogram(minutes(timeMT-time));
xlabel('Minutes');
ylabel('PDF');

sigmaTI = 20;
timeTI = time + minutes(random('Weibull',sigmaTI,2,length(time),1));

X4MT = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeMT));
eMT = X4MT-X4;

X4TI = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(timeTI));
eTI = X4TI-X4; 

%% Error distribution X* = X + XV; what is V? 

eNorm =(eMT)./X4;


[feNgXS,XeNS,YeNS] = conditional_pdf(eNorm,X4,-20:0.5:20,0:1:20);
figure;
p=pcolor(XeNS,YeNS,feNgXS);
set(p,'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X_{geotail}');
xlabel('P( X^*-X/X|X_{geotail})');
title('Normalized Error');
colormap(inferno);
xlim([-1,1]);

figure; histogram(eNorm(~isnan(eNorm)),-5:0.1:5,'Normalization','pdf');

%%
pdf1 = fitdist(eNorm(~isnan(eNorm) & X4>2 & X4<30),'kernel');
figure;
x = -5:0.1:5;
histogram(eNorm(~isnan(eNorm)),x,'Normalization','pdf');
hold on;
histogram(random(pdf1,1000,1),x,'Normalization','pdf');
hold on;
plot(x,pdf(pdf1,x));
hold on;
histogram(eNorm(~isnan(eNorm) & X4>3 & X4<20),x,'Normalization','pdf');

%%
% Xm_3c = X + X.*random(pdf1,Nsamples,Nensemble);
Xm_3b = X + X.*random(pdf1,Nsamples,Nensemble);
eNormT = (Xm_3b - X)./X;

[feNgXT,XeNT,YeNT] = conditional_pdf(eNormT,X,-20:0.5:20,0:1:20);


figure;
p=pcolor(XeNT,YeNT,feNgXT);
set(p,'EdgeColor','none');
colorbar;
caxis([0,1]);
ylabel('X');
xlabel('P( (X^*-X)/X|X)');
title('Normalized Error');
colormap(inferno);
xlim([-1,1]);

figure; histogram(eNormT,-5:0.1:5,'Normalization','pdf');




%%
filterS = X4<6 & X4>5;
filterT = X<6 & X>5;

figure; 
histogram(e1(filterT),-5:0.1:5,'Normalization','pdf'); hold on;
histogram(e(filterT),-5:0.1:5,'Normalization','pdf'); hold on; 
histogram(eS(filterS),-5:0.1:5,'Normalization','pdf'); xlim([-5,+5]);
legend('Norm','Lognorm','Data');

% ylim([0,1]);
    

%%
figure; 
plot(YMTS(:,1)',(std(feMTgXS')./YMTS(:,1)')); ylim([0,0.1]);
hold on;
plot(YTIS(:,1)',(std(feTIgXS')./YTIS(:,1)')); ylim([0,0.1]);
hold on;
plot(YTS(:,1)',(std(feTgXS')./YTS(:,1)')); ylim([0,0.1]);

%%
i=400; 
f1 = fit(XTS(i,:).',feTgXS(i,:).','gauss1'); 
figure; plot(XTS(i,:),(feTgXS(i,:))); 
hold on; 
plot(XTS(i,:),f1(XTS(i,:))); 
title(['X = ',num2str(YTS(i,1)),'; \sigma_{TI} = ' num2str(f1.c1),'; Ratio =',num2str(f1.c1./YTS(i,1))]);

f2 = fit(XMTS(i,:).',feMTgXS(i,:).','gauss1'); 
figure; plot(XMTS(i,:),(feMTgXS(i,:))); 
hold on; 
plot(XMTS(i,:),f2(XMTS(i,:))); 
title(['X = ',num2str(YMTS(i,1)),'; \sigma_{TM} = ' num2str(f2.c1),'; Ratio =',num2str(f2.c1./YMTS(i,1))]);

%% Functions

function curve = create_curve(X, Y, Ei)

if nargin<3
    Ei = 100;
end
Y = Y(:);
Y(Y==999999)=nan;

X1 = X(~isnan(X) & ~isnan(Y));
Y1 = Y(~isnan(X) & ~isnan(Y));
X = X1;
Y = Y1;

[xindx, E] = discretize(X(:),Ei);

for i = 1:max(xindx)
    curve.YgX(i) = nanmean(Y(xindx==i));
    curve.stdYgX(i) = nanstd(Y(xindx==i));
    curve.NSamples(i) = sum(xindx==i & ~isnan(Y));
    curve.SEM(i) = nanstd(Y(xindx==i))./sqrt(curve.NSamples(i));
    curve.ts(i,:) = tinv([0.025 0.975],curve.NSamples(i)-1);
    curve.CI(i,:) = curve.YgX(i) + curve.ts(i,:)*curve.SEM(i);
    curve.XBins(i) = 0.5*(E(i)+E(i+1));
end

curve.E = E; 

end

function [fXgY,XX,YY] = conditional_pdf(X, Y, gridx, gridy)

X1 = X(~isnan(X) & ~isnan(Y));
Y1 = Y(~isnan(X) & ~isnan(Y));
X = X1;
Y = Y1;

sz = 2^10;
Y = Y(:);
Y(Y==999999)=nan;

[bandwidth,fXY,XX,YY]=kde2d([X,Y],sz,[min(gridx),min(gridy)],[max(gridx),max(gridy)]);
[fY, YY1] = ksdensity(Y,YY(:,1));

fXgY = fXY./repmat(fY,1,sz);

% [XX, YY] = meshgrid(gridx,gridy);
% 
% [fxy, yii] = ksdensity([X,Y],[XX(:),YY(:)]);
% Fxy = scatteredInterpolant(yii(:,1),yii(:,2),fxy);
% 
% [fy, yi] = ksdensity(Y,gridy);
% Fy = griddedInterpolant(yi,fy);
% 
% cpdf = Fxy(XX,YY)./(Fy(gridy)')./(sum((Fxy(XX,YY)./(Fy(gridy)')),2));

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

function [bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
% fast and accurate state-of-the-art
% bivariate kernel density estimator
% with diagonal bandwidth matrix.
% The kernel is assumed to be Gaussian.
% The two bandwidth parameters are
% chosen optimally without ever
% using/assuming a parametric model for the data or any "rules of thumb".
% Unlike many other procedures, this one
% is immune to accuracy failures in the estimation of
% multimodal densities with widely separated modes (see examples).
% INPUTS: data - an N by 2 array with continuous data
%            n - size of the n by n grid over which the density is computed
%                n has to be a power of 2, otherwise n=2^ceil(log2(n));
%                the default value is 2^8;
% MIN_XY,MAX_XY- limits of the bounding box over which the density is computed;
%                the format is:
%                MIN_XY=[lower_Xlim,lower_Ylim]
%                MAX_XY=[upper_Xlim,upper_Ylim].
%                The dafault limits are computed as:
%                MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
%                MAX_XY=MAX+Range/4; MIN_XY=MIN-Range/4;
% OUTPUT: bandwidth - a row vector with the two optimal
%                     bandwidths for a bivaroate Gaussian kernel;
%                     the format is:
%                     bandwidth=[bandwidth_X, bandwidth_Y];
%          density  - an n by n matrix containing the density values over the n by n grid;
%                     density is not computed unless the function is asked for such an output;
%              X,Y  - the meshgrid over which the variable "density" has been computed;
%                     the intended usage is as follows:
%                     surf(X,Y,density)
% Example (simple Gaussian mixture)
% clear all
%   % generate a Gaussian mixture with distant modes
%   data=[randn(500,2);
%       randn(500,1)+3.5, randn(500,1);];
%   % call the routine
%     [bandwidth,density,X,Y]=kde2d(data);
%   % plot the data and the density estimate
%     contour3(X,Y,density,50), hold on
%     plot(data(:,1),data(:,2),'r.','MarkerSize',5)
%
% Example (Gaussian mixture with distant modes):
%
% clear all
%  % generate a Gaussian mixture with distant modes
%  data=[randn(100,1), randn(100,1)/4;
%      randn(100,1)+18, randn(100,1);
%      randn(100,1)+15, randn(100,1)/2-18;];
%  % call the routine
%    [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none'), view([0,60])
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
%
% Example (Sinusoidal density):
%
% clear all
%   X=rand(1000,1); Y=sin(X*10*pi)+randn(size(X))/3; data=[X,Y];
%  % apply routine
%  [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none'), view([0,70])
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
%
%  Reference:
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
global N A2 I
if nargin<2
    n=2^8;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;
N=size(data,1);
if nargin<3
    MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
    MAX_XY=MAX+Range/2; MIN_XY=MIN-Range/2;
end
scaling=MAX_XY-MIN_XY;
if N<=size(data,2)
    error('data has to be an N by 2 array where each row represents a two dimensional observation')
end
transformed_data=(data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
%bin the data uniformly using regular grid;
initial_data=ndhist(transformed_data,n);
% discrete cosine transform of initial data
a= dct2d(initial_data);
% now compute the optimal bandwidth^2
  I=(0:n-1).^2; A2=a.^2;
 t_star=root(@(t)(t-evolve(t)),N);
p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
% smooth the discrete cosine transform of initial data using t_star
a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
% now apply the inverse discrete cosine transform
if nargout>1
    density=idct2d(a_t)*(numel(a_t)/prod(scaling));
	density(density<0)=eps; % remove any negative density values
    [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
end
bandwidth=sqrt([t_x,t_y]).*scaling; 
end
%#######################################
function  [out,time]=evolve(t)
global N
Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
time=(2*pi*N*Sum_func)^(-1/3);
out=(t-time)/time;
end
%#######################################
function out=func(s,t)
global N
if sum(s)<=4
    Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
    time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
    out=psi(s,time);
else
    out=psi(s,t);
end
end
%#######################################
function out=psi(s,Time)
global I A2
% s is a vector
w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
wx=w.*(I.^s(1));
wy=w.*(I.^s(2));
out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
end
%#######################################
function out=K(s)
out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
end
%#######################################
function data=dct2d(data)
% computes the 2 dimensional discrete cosine transform of data
% data is an nd cube
[nrows,ncols]= size(data);
if nrows~=ncols
    error('data is not a square array!')
end
% Compute weights to multiply DFT coefficients
w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
weight=w(:,ones(1,ncols));
data=dct1d(dct1d(data)')';
    function transform1d=dct1d(x)
        % Re-order the elements of the columns of x
        x = [ x(1:2:end,:); x(end:-2:2,:) ];
        % Multiply FFT by weights:
        transform1d = real(weight.* fft(x));
    end
end
%#######################################
function data = idct2d(data)
% computes the 2 dimensional inverse discrete cosine transform
[nrows,ncols]=size(data);
% Compute wieghts
w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
weights=w(:,ones(1,ncols));
data=idct1d(idct1d(data)');
    function out=idct1d(x)
        y = real(ifft(weights.*x));
        out = zeros(nrows,ncols);
        out(1:2:nrows,:) = y(1:nrows/2,:);
        out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
    end
end
%#######################################
function binned_data=ndhist(data,M)
% this function computes the histogram
% of an n-dimensional data set;
% 'data' is nrows by n columns
% M is the number of bins used in each dimension
% so that 'binned_data' is a hypercube with
% size length equal to M;
[nrows,ncols]=size(data);
bins=zeros(nrows,ncols);
for i=1:ncols
    [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
    bins(:,i) = min(bins(:,i),M);
end
% Combine the  vectors of 1D bin counts into a grid of nD bin
% counts.
binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t=root(f,N)
% try to find smallest root whenever there is more than one
N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
tol=10^-12+0.01*(N-50)/1000;
flag=0;
while flag==0
    try
        t=fzero(f,[0,tol]);
        flag=1;
    catch
        tol=min(tol*2,.1); % double search interval
    end
    if tol==.1 % if all else fails
        t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
    end
end
end

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

