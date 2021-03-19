%% Routine that constructs a super-posed epoch analysis of substorms
% All substorms in the database, and also ones that are observed at PFISR

% Date: 10th February 2021

% Initialization

if exist('prevSuperMagFileType','var')
    if prevSuperMagFileType ~= superMagFileType
        clear all;
    end
end

storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\'; 
superMagFileStr = [storeDir,'substorms_superMag_20201130.txt'];
superMagFileType = 1; % 1 -> superMag database, 2 -> Colin Forsyth's data
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
workDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir';

% Time window around onset to process and plot data
preOnsetWindow = duration(3,0,0); 
postOnsetWindow = duration(3,0,0);

if ~exist('prevSuperMagFileType','var')
    prevSuperMagFileType = 0;
end

% Constructing table of substorms while tracking PFISR experiments
if ~(exist('T','var') && prevSuperMagFileType == superMagFileType) 
    % runs only when the program starts for the first time
    timeMinStr = "01 Dec 2006";
    timeMaxStr = "31 Dec 2019";
    [T, superMag] = substorm_create_table(superMagFileStr, superMagFileType,...
      timeMinStr, timeMaxStr);
    
else
    warning('Substorm table already exists. See table T.');
end

%%

if prevSuperMagFileType == 0 % For the first run of this routine
    
    % Cataloging stored substorms
    fileNameList = struct2cell(dir([workDir,strcat(filesep,'*_pfisrData.h5')]));
    filePathStr = strcat(strcat(workDir,filesep),string(fileNameList(1,:)'));
    for i=1:1:length(filePathStr)
        storedStormID(:,i) = str2double(fileNameList{1,i}(1:regexp(fileNameList{1,i},'_')-1));
    end
    
    % Loading all relevant data variables
    omni = extract_omni_data(omniFile);
    

else  
    warning(['Since I think this is not the first run,',...
        ' did not load data variables again.',...
        ' Set prevSuperMagFileType=0 to rerun this section.']);
end


%% Identifying and defining isolated substorms

% Identify peak SML and peak SML Time during the substorm window 
% Substorm windows - 1 hr and + 30 minutes from the onset
% (a proxy for the strength of the substorm)

if ~(exist('T','var') && prevSuperMagFileType == superMagFileType) 
    
    for i=1:1:length(T.stormID)
     
        if sum(storedStormID==T.stormID(i))>0
            T.storageLocation(i) = filePathStr(storedStormID==T.stormID(i));
        end
    
    tempTime = T.Time(i)-duration(1,0,0):duration(0,1,0):T.Time(i)+duration(0,30,0);
    [T.peakSML(i),indx2] = min(omni.Fsml(datenum(tempTime)));
    T.peakSMLTime(i)= tempTime(indx2);
    
    
    
    end

    T = identify_qualified_substorm_duration(T); 
    
else
    
    warning('Substorm table already exists. See table T.');

end

%% Loading OMNI data corresponding to each substorm onto the table

if ~(exist('T','var') && prevSuperMagFileType == superMagFileType) 
    T = add_omni_array_to_table(T, omni, preOnsetWindow, postOnsetWindow);
else   
    warning('Substorm table already exists, so not updating data. See table T.');
end

%% Loading PFISR data corresponding to each substorm with pfisr data
tic
if ~(exist('T','var') && prevSuperMagFileType == superMagFileType) 
    T=add_pfisr_array_to_table(T);
else
    warning('Substorm table already exists, so not updating data. See table T.');
end
toc

%% Plotting

condition = T.previousSubstormDuration>duration(3,0,0)... %    
    & T.nextSubstormDuration>duration(3,0,0);%...
%     & absDiffMLT(T.MLT, T.PFISR_MLT)<2 ...
%     & T.MLAT<66;
%     & ~ismissing(T.storageLocation); %    


plot_superposed_epoch(T, condition, {'BzArray','EArray','symHArray','pfisrNe'});
% plot_epoch(T, 30581, {'BzArray','EArray','symHArray','pfisrHall','pfisrNe','pfisrEflux'});


%% Rough work

% Distribution of these substorm classes with respect to amount of energy
% accumulated and released 
multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0) & T.peakSML<-500 & T.peakSML>-600;
singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0) & T.peakSML<-500 & T.peakSML>-600;


figure; 
edges = logspace(3,6,100);
subplot(2,1,1)
histogram(trapz(minutes(T.timeArrayRelOnset(1,1:180)),interp_nans(T.EArray(multiSubstorm,1:180)')',2),edges,'Normalization','countdensity'); 
hold on; 
histogram(trapz(minutes(T.timeArrayRelOnset(1,1:180)),interp_nans(T.EArray(singleSubstorm,1:180)')',2),edges,'Normalization','countdensity'); 
set(gca,'XScale','log');
xlabel('\int (VxB)_{sw} before Onset'); 
legend('Multi-onset','Single-onset');
xlim([10^3, 10^6]);

subplot(2,1,2)
histogram(trapz(minutes(T.timeArrayRelOnset(1,180:361)),interp_nans(T.EArray(multiSubstorm,180:361)')',2),edges,'Normalization','countdensity'); 
hold on; 
histogram(trapz(minutes(T.timeArrayRelOnset(1,180:361)),interp_nans(T.EArray(singleSubstorm,180:361)')',2),edges,'Normalization','countdensity'); 
set(gca,'XScale','log');
xlabel('\int (VxB)_{sw} after Onset'); 
legend('Multi-onset','Single-onset');
xlim([10^3, 10^6]);

%% Averaging a set of substorms before plotting the distribution 
EArray1 = T.EArray; 
% EArray0(EArray0<-1) = 0; 
% EArray1 = interp_nans(T.EArray')';
%%
clearvars multiStorm singleStorm neitherStorm
SML_max = -100;
SML_min = -900;
multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0) & T.peakSML<SML_max & T.peakSML>SML_min;
singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0) & T.peakSML<SML_max & T.peakSML>SML_min;
% multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0);
% singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0);
neitherSubstorm = ~(multiSubstorm | singleSubstorm) & T.peakSML<SML_max & T.peakSML>SML_min;

multiStormID = T.stormID(multiSubstorm);
singleStormID = T.stormID(singleSubstorm);
neitherStormID = T.stormID(neitherSubstorm);
sampleSize = 100;
ensembleSize = 5000;
for i= 1:ensembleSize
    multiStorm.EnsembleID(i,:)=multiStormID(randperm(length(multiStormID),sampleSize));
    singleStorm.EnsembleID(i,:)=singleStormID(randperm(length(singleStormID),sampleSize));
    neitherStorm.EnsembleID(i,:)=neitherStormID(randperm(length(neitherStormID),sampleSize));
    
    multiStorm.EArray(i,:) = nanmean((EArray1(multiStorm.EnsembleID(i,:),:)),1);
    singleStorm.EArray(i,:) = nanmean((EArray1(singleStorm.EnsembleID(i,:),:)),1);
    neitherStorm.EArray(i,:) = nanmean((EArray1(neitherStorm.EnsembleID(i,:),:)),1);
    
end
singleStorm.timeArrayRelOnset = T.timeArrayRelOnset(1,:);

growthTIndex = 60:180; %Minutes
expansionTIndex = 181:301; %Minutes
%%
figure; 
edges = logspace(3,8,200);
subplot(2,1,1)

h1 = histogram(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),singleStorm.EArray(:,growthTIndex),2),edges,'Normalization','countdensity'); 
hold on;
h2 = histogram(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),multiStorm.EArray(:,growthTIndex),2),edges,'Normalization','countdensity'); 
hold on;
h3 = histogram(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),neitherStorm.EArray(:,growthTIndex),2),edges,'Normalization','countdensity'); 
 
set(gca,'XScale','log');
xlabel('$\int_{t_0-2 \mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
title('Pre-Onset');
ylabel('Count Density');
legend(['Single-onset # ',num2str(length(singleStormID))],...
    ['Multi-onset # ',num2str(length(multiStormID))],...
    ['Non-isolated # ',num2str(length(neitherStormID))]);
xlim([10^3, 10^8]);
text(0.01,0.3,{['#Substorms-per-Ensemble: ',num2str(sampleSize)],['#Ensembles: ',num2str(ensembleSize)]},'Units','Normalized');
text(0.01,0.9,{[num2str(SML_max),'nT > SML > ',num2str(SML_min),'nT']},'Units','Normalize');
% ylim([0,0.15]);


subplot(2,1,2)
histogram(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),singleStorm.EArray(:,expansionTIndex),2),edges,'Normalization','countdensity'); 
hold on;
histogram(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),multiStorm.EArray(:,expansionTIndex),2),edges,'Normalization','countdensity'); 
hold on;
histogram(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),neitherStorm.EArray(:,expansionTIndex),2),edges,'Normalization','countdensity'); 
 

set(gca,'XScale','log');
xlabel('$\int_{t_0}^{t_0+2\mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 

legend(['Single-onset # ',num2str(length(singleStormID))],...
    ['Multi-onset # ',num2str(length(multiStormID))],...
    ['Non-isolated # ',num2str(length(neitherStormID))]);

title('Post-Onset');
ylabel('Count Density');
xlim([10^3, 10^8]);
% ylim([0,0.15]);
%E = 1e-9.*1e7.*(velocity.*10.^3).*((B.*10^-9).^2).*(sin(theta_c/2)).^4*l_0^2; %GW 

figure;
tot = h1.Values + h2.Values + h3.Values;
plot(h1.BinEdges(1:end-1) , h1.Values./tot);
hold on;
plot(h2.BinEdges(1:end-1) , h2.Values./tot);
hold on;
plot(h3.BinEdges(1:end-1) , h3.Values./tot);
set(gca,'XScale','log');
ylabel('Probability');
xlabel('$\int_{t_0-2 \mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 

%% Plot of correlation between pre-onset and post-onset driving
clearvars multiStorm singleStorm neitherStorm
cmap = inferno;
growthTIndex = (60:180); %Minutes
expansionTIndex = (181:361); %Minutes
% EArray2 = interp_nans(T.EArray')';
EArray2 = T.EArray;

SML = -100:-100:-1200;

j = 3;

SML_max = SML(j);
SML_min = SML(j+1);
preOnsetE = nansum(T.EArray(:,growthTIndex),2);
postOnsetE = nansum(T.EArray(:,expansionTIndex),2);
postOMin = 10^4;
postOMax = 5*10^4;
% postOMin = 10^2;
% postOMax = 10^6;

multiSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
 T.nextSubstormDuration<duration(3,0,0)...
 & T.peakSML<SML_max & T.peakSML>SML_min & sum(~isnan(T.EArray),2)>320&...
 postOnsetE<postOMax & postOnsetE>postOMin;

singleSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
    T.nextSubstormDuration>duration(3,0,0) & T.peakSML<SML_max &...
    T.peakSML>SML_min & sum(~isnan(T.EArray),2)>320 &...
    postOnsetE<postOMax & postOnsetE>postOMin;

edges = logspace(3,5,50);

figure; 
x1 = nansum(EArray2(singleSubstorm,growthTIndex),2);
y1 = nansum(EArray2(singleSubstorm,expansionTIndex),2);
scatter(x1,y1,'.','MarkerEdgeColor',cmap(1,:));
b1 = [ones(length(x1),1) x1]\y1;
x = logspace(2,6,100);
hold on;
plot(x,x*b1(2)+b1(1),'Color',cmap(1,:));
c1 = corrcoef(x1,y1);
text(0.7,0.6,['r_0_0 = ', num2str(c1(1,2),2)],'Units','Normalized','Color',cmap(1,:));
set(gca,'YScale','log','XScale','log','XLim',[10^2,10^6],'YLim',[10^2,10^6]);
xlabel('Pre-onset: $\int_{t_0-2\mathrm{Hr}}^{t_0 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
ylabel('Post-onset: $\int_{t_0}^{t_0+2 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
hold on; 

x2 = nansum(EArray2(multiSubstorm,growthTIndex),2);
y2 = nansum(EArray2(multiSubstorm,expansionTIndex),2);
scatter(x2,y2,'.','MarkerEdgeColor',cmap(160,:));
b2 = [ones(length(x2),1) x2]\y2;
hold on;
plot(x,x*b2(2)+b2(1),'Color',cmap(160,:));
c2 = corrcoef(x2,y2);
text(0.7,0.9,['r_0_0 = ', num2str(c2(1,2),2)],'Units','Normalized','Color',cmap(160,:));

legend('Single-onset','Single-onset linear fit','Multi-onset','Multi-onset linear fit','Location','southeast');
text(0.1,0.9,[num2str(SML(j)),'<SML<',num2str(SML(j+1)),'nT'],'Units','Normalized');

figure;
histogram(x1,edges,'Normalization','countdensity','FaceColor',cmap(1,:)); hold on;
histogram(x2,edges,'Normalization','countdensity','FaceColor',cmap(160,:));
set(gca,'XScale','log');
xlabel('Pre-onset: $\int_{t_0-2\mathrm{Hr}}^{t_0 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
ylabel('Count density');
text(0.1,0.9,[num2str(SML(j)),'<SML<',num2str(SML(j+1)),'nT'],'Units','Normalized');
legend('Single-onset','Multi-onset');

figure;
histogram(bootstrp(5000,@mean,x1),'Normalization','countdensity','FaceColor',cmap(1,:)); hold on;
histogram(bootstrp(5000,@mean,x2),'Normalization','countdensity','FaceColor',cmap(160,:));
set(gca,'XScale','log');
xlabel('Pre-onset: $ \langle \int_{t_0-2\mathrm{Hr}}^{t_0 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t \rangle_{substorm-samples}$ [GJ]','Interpreter','latex'); 
ylabel('Count density');
text(0.1,0.9,[num2str(SML(j)),'<SML<',num2str(SML(j+1)),'nT'],'Units','Normalized');
legend('Single-onset','Multi-onset');

%% Plot of Driving vs SML (without averaging)
clearvars multiStorm singleStorm neitherStorm
cmap = inferno;
% growthTIndex = (60:180); %Minutes
growthTIndex = 12;
expansionTIndex = (181:361); %Minutes

xAxisRange = [1,10^4];
yAxisRange = [0,1400];
% xLabelStr = 'Pre-onset: $\int_{t_0-2\mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]';
xLabelStr = ['$V B^2 \sin^4 \frac{\theta_c}{2} L_0^2$ [GJ] at $t=t_0+$',num2str(growthTIndex-180),'min ']; 
% yLabelStr = 'Peak $|SML|$ [nT]';
yLabelStr = '$|SML|$ [nT]';

xAxisBootRange = [1,1000];
% xLabelBootStr = 'Pre-onset: $\langle \int_{t_0-2\mathrm{Hr}}^{t_0 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t \rangle_{substorm-samples}$ [GJ]';
xLabelBootStr = ['$\langle V B^2 \sin^4 \frac{\theta_c}{2} L_0^2\rangle_{substorm-samples}$ [GJ] at $t=t_0+$',num2str(growthTIndex-180),'min '];
% yLabelBootStr = '$\langle \mathrm{Peak \ |SML|} \rangle_{substorm-samples}$ [nT]';
yLabelBootStr = yLabelStr;

% EArray2 = interp_nans(T.EArray')';
EArray2 = T.EArray;
peakSML = T.smlArray(:,growthTIndex); 
% preOnsetE = nansum(T.EArray(:,growthTIndex),2);
preOnsetE = nansum(T.EArray(:,growthTIndex),2);
postOnsetE = nansum(T.EArray(:,expansionTIndex),2);

multiSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
 T.nextSubstormDuration<duration(3,0,0)...
 & sum(~isnan(T.EArray),2)>320;

singleSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
    T.nextSubstormDuration>duration(3,0,0) & sum(~isnan(T.EArray),2)>320;

SML = -100:-50:-1200;

for j=1:length(SML)-1
    
    SML_max = SML(j);
    SML_min = SML(j+1);
    multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0) & peakSML<SML_max & peakSML>SML_min;
    singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0) & peakSML<SML_max & peakSML>SML_min;
    % multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0);
    % singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0);
    neitherSubstorm = ~(multiSubstorm | singleSubstorm) & peakSML<SML_max & peakSML>SML_min;

    multiStorm.index(j,:) = multiSubstorm;
    singleStorm.index(j,:) = singleSubstorm;
    neitherStorm.index(j,:) = neitherSubstorm;
    
end

[~,singleStormID]=(ind2sub(size(singleStorm.index),find(singleStorm.index)));
[~,multiStormID]=(ind2sub(size(multiStorm.index),find(multiStorm.index)));
x1 = preOnsetE(singleStormID);
% x11 = EArray2(singleSubstorm,180);
y1 = peakSML(singleStormID);
x2 = preOnsetE(multiStormID);
% x22 = EArray2(multiSubstorm,180);
y2 = peakSML(multiStormID);

figure; 
scatter(x1,abs(y1),10,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.4);
hold on;
scatter(x2,abs(y2),10,'filled','MarkerFaceColor',cmap(160,:),'MarkerFaceAlpha',0.4);
% b1 = [ones(length(x1),1) x1]\y1;
x = logspace(2,6,100);
% hold on;
% plot(x,x*b1(2)+b1(1),'Color',cmap(1,:));
% c1 = corrcoef(x1,y1);
% text(0.7,0.6,['r_0_0 = ', num2str(c1(1,2),2)],'Units','Normalized','Color',cmap(1,:));
set(gca,'YScale','linear','XScale','log','XLim',xAxisRange,'YLim',yAxisRange);
xlabel(xLabelStr,'Interpreter','latex'); 
ylabel(yLabelStr,'Interpreter','latex'); 
legend('Single-onset','Multi-onset');

% Bootstrap the above
nBoot = 1000;
for iSML = 1:length(SML)-1
    try
        singleStorm.preOnsetE(iSML,:) = bootstrp(nBoot,@mean,preOnsetE(singleStorm.index(iSML,:)));
        singleStorm.postOnsetE(iSML,:) = bootstrp(nBoot,@mean,postOnsetE(singleStorm.index(iSML,:)));
    catch
        singleStorm.preOnsetE(iSML,:)=nan(1,nBoot);
        singleStorm.postOnsetE(iSML,:)=nan(1,nBoot);
    end
    
    try
        singleStorm.peakSML(iSML,:) = bootstrp(nBoot,@mean,peakSML(singleStorm.index(iSML,:)));
    catch
        singleStorm.peakSML(iSML,:) = nan(1,nBoot);
    end
    
    try
        multiStorm.preOnsetE(iSML,:) = bootstrp(nBoot,@mean,preOnsetE(multiStorm.index(iSML,:)));
        multiStorm.postOnsetE(iSML,:) = bootstrp(nBoot,@mean,postOnsetE(multiStorm.index(iSML,:)));
    catch
        multiStorm.preOnsetE(iSML,:) = nan(1,nBoot);
        multiStorm.postOnsetE(iSML,:) = nan(1,nBoot);
    end
    
    try
        multiStorm.peakSML(iSML,:) = bootstrp(nBoot,@mean,peakSML(multiStorm.index(iSML,:))); 
    catch
        multiStorm.peakSML(iSML,:) = nan(1,nBoot);
    end
    
end

figure; 
x1 = singleStorm.preOnsetE(:);
xm1 = mean(singleStorm.preOnsetE,2);

y1 = singleStorm.peakSML(:);
ym1 = mean(singleStorm.peakSML,2);

x2 = multiStorm.preOnsetE(:);
xm2 = mean(multiStorm.preOnsetE,2);
y2 = multiStorm.peakSML(:);
ym2 = mean(multiStorm.peakSML,2);

p1=scatter(x1,abs(y1),10,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.03);
hold on;

p2=scatter(x2,abs(y2),10,'filled','MarkerFaceColor',cmap(160,:),'MarkerFaceAlpha',0.03);
hold on;
plot(xm1,abs(ym1),'Color',cmap(60,:),'LineWidth',2);
hold on;
plot(xm2,abs(ym2),'Color',cmap(200,:),'LineWidth',2);
x = logspace(2,6,100);
set(gca,'YScale','linear','XScale','linear','XLim',xAxisBootRange,'YLim',yAxisRange);
xlabel(xLabelBootStr,'Interpreter','latex'); 
ylabel(yLabelBootStr,'Interpreter','latex'); 
legend([p1,p2],{'Single-onset','Multi-onset'});
%%
% Post-onset
figure; 
x1 = singleStorm.postOnsetE(:);
xm1 = mean(singleStorm.postOnsetE,2);

y1 = singleStorm.peakSML(:);
ym1 = mean(singleStorm.peakSML,2);

x2 = multiStorm.postOnsetE(:);
xm2 = mean(multiStorm.postOnsetE,2);
y2 = multiStorm.peakSML(:);
ym2 = mean(multiStorm.peakSML,2);

p1=scatter(x1,y1,'.','MarkerEdgeColor',cmap(1,:));
hold on;

p2=scatter(x2,y2,'.','MarkerEdgeColor',cmap(160,:));
hold on;
plot(xm1,ym1,'Color',cmap(60,:),'LineWidth',2);
hold on;
plot(xm2,ym2,'Color',cmap(200,:),'LineWidth',2);
x = logspace(2,6,100);
set(gca,'YScale','linear','XScale','log','XLim',[5*10^3,10^6],'YLim',[-1400,0]);
xlabel('Post-onset: $\langle \int_{t_0}^{t_0 +2\mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t \rangle_{substorm-samples}$ [GJ]','Interpreter','latex'); 
ylabel('$\langle \mathrm{Peak \ SML} \rangle_{substorm-samples}$ [nT]','Interpreter','latex'); 
legend([p1,p2],{'Single-onset','Multi-onset'});


%% Plot of Driving (EMF) vs SML (current_estimate) to demonstrate linearity 
clearvars multiStorm singleStorm neitherStorm
cmap = inferno;
% growthTIndex = (60:180); %Minutes
growthTIndex = 300;
expansionTIndex = (181:361); %Minutes

xAxisRange = [0,10];
yAxisRange = [0,1400];
% xLabelStr = 'Pre-onset: $\int_{t_0-2\mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]';
xLabelStr = ['$E_M(V,B_T,\theta_c)$ [mV/m] at $t=t_0+$',num2str(growthTIndex-180),'min ']; 
% yLabelStr = 'Peak $|SML|$ [nT]';
yLabelStr = '$J_{iono}(SML)$ [mA/m$^2$]';

xAxisBootRange = [0,5];
% xLabelBootStr = 'Pre-onset: $\langle \int_{t_0-2\mathrm{Hr}}^{t_0 \mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t \rangle_{substorm-samples}$ [GJ]';
xLabelBootStr = ['$\langle E_M \rangle_{substorm-samples}$ [mV/m] at $t=t_0+$',num2str(growthTIndex-180),'min '];
% yLabelBootStr = '$\langle \mathrm{Peak \ |SML|} \rangle_{substorm-samples}$ [nT]';
yLabelBootStr = yLabelStr;

% EArray2 = interp_nans(T.EArray')';
EArray2 = T.EArray;
C=define_universal_constants;
peakSML = (2./C.mu0)*T.smlArray(:,growthTIndex).*10^-6;%[mA/m^2 (?)]

% preOnsetE = nansum(T.EArray(:,growthTIndex),2);
preOnsetE = nansum(T.eklArray(:,growthTIndex),2).*10^-3; %[mV/m] See Kan and 
postOnsetE = nansum(T.eklArray(:,expansionTIndex),2).*10^-3;

multiSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
 T.nextSubstormDuration<duration(3,0,0)...
 & sum(~isnan(T.EArray),2)>320;

singleSubstorm = T.previousSubstormDuration>duration(3,0,0) &...
    T.nextSubstormDuration>duration(3,0,0) & sum(~isnan(T.EArray),2)>320;

SML = -100:-50:-1200;

for j=1:length(SML)-1
    
    SML_max = SML(j);
    SML_min = SML(j+1);
    multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0) & peakSML<SML_max & peakSML>SML_min;
    singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0) & peakSML<SML_max & peakSML>SML_min;
    % multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0);
    % singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0);
    neitherSubstorm = ~(multiSubstorm | singleSubstorm) & peakSML<SML_max & peakSML>SML_min;

    multiStorm.index(j,:) = multiSubstorm;
    singleStorm.index(j,:) = singleSubstorm;
    neitherStorm.index(j,:) = neitherSubstorm;
    
end

[~,singleStormID]=(ind2sub(size(singleStorm.index),find(singleStorm.index)));
[~,multiStormID]=(ind2sub(size(multiStorm.index),find(multiStorm.index)));
x1 = preOnsetE(singleStormID);
% x11 = EArray2(singleSubstorm,180);
y1 = peakSML(singleStormID);
x2 = preOnsetE(multiStormID);
% x22 = EArray2(multiSubstorm,180);
y2 = peakSML(multiStormID);

figure; 
scatter(x1,abs(y1),10,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.4);
hold on;
scatter(x2,abs(y2),10,'filled','MarkerFaceColor',cmap(160,:),'MarkerFaceAlpha',0.4);
% b1 = [ones(length(x1),1) x1]\y1;
x = logspace(2,6,100);
% hold on;
% plot(x,x*b1(2)+b1(1),'Color',cmap(1,:));
% c1 = corrcoef(x1,y1);
% text(0.7,0.6,['r_0_0 = ', num2str(c1(1,2),2)],'Units','Normalized','Color',cmap(1,:));
set(gca,'YScale','linear','XScale','linear','XLim',xAxisRange,'YLim',yAxisRange);
xlabel(xLabelStr,'Interpreter','latex'); 
ylabel(yLabelStr,'Interpreter','latex'); 
legend('Single-onset','Multi-onset');

% Bootstrap the above
nBoot = 5000;
for iSML = 1:length(SML)-1
    try
        singleStorm.preOnsetE(iSML,:) = bootstrp(nBoot,@mean,preOnsetE(singleStorm.index(iSML,:)));
        singleStorm.postOnsetE(iSML,:) = bootstrp(nBoot,@mean,postOnsetE(singleStorm.index(iSML,:)));
    catch
        singleStorm.preOnsetE(iSML,:)=nan(1,nBoot);
        singleStorm.postOnsetE(iSML,:)=nan(1,nBoot);
    end
    
    try
        singleStorm.peakSML(iSML,:) = bootstrp(nBoot,@mean,peakSML(singleStorm.index(iSML,:)));
    catch
        singleStorm.peakSML(iSML,:) = nan(1,nBoot);
    end
    
    try
        multiStorm.preOnsetE(iSML,:) = bootstrp(nBoot,@mean,preOnsetE(multiStorm.index(iSML,:)));
        multiStorm.postOnsetE(iSML,:) = bootstrp(nBoot,@mean,postOnsetE(multiStorm.index(iSML,:)));
    catch
        multiStorm.preOnsetE(iSML,:) = nan(1,nBoot);
        multiStorm.postOnsetE(iSML,:) = nan(1,nBoot);
    end
    
    try
        multiStorm.peakSML(iSML,:) = bootstrp(nBoot,@mean,peakSML(multiStorm.index(iSML,:))); 
    catch
        multiStorm.peakSML(iSML,:) = nan(1,nBoot);
    end
    
end
%%
figure; 
x=linspace(xAxisRange(1),xAxisRange(2),10);
x1 = singleStorm.preOnsetE(:);
xm1 = mean(singleStorm.preOnsetE,2);

y1 = abs(singleStorm.peakSML(:));
ym1 = abs(mean(singleStorm.peakSML,2));



x2 = multiStorm.preOnsetE(:);
xm2 = mean(multiStorm.preOnsetE,2);
y2 = abs(multiStorm.peakSML(:));
ym2 = abs(mean(multiStorm.peakSML,2));

% fit
[f1,gof1] = fit(x1(~isnan(x1)),y1(~isnan(x1)),'poly1','Robust','Bisquare','Exclude', y1(~isnan(x1)) > 600);
[f2,gof2] = fit(x2(~isnan(x2)),y2(~isnan(x2)),'poly1','Robust','Bisquare','Exclude', y2(~isnan(x2)) > 1200);

% plot
scatter(x1,y1,10,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.03);
hold on;

scatter(x2,y2,10,'filled','MarkerFaceColor',cmap(160,:),'MarkerFaceAlpha',0.03);
hold on;
p1=plot(x,f1(x),'Color',cmap(60,:),'LineWidth',2);
hold on;
p2=plot(x,f2(x),'Color',cmap(200,:),'LineWidth',2);
x = logspace(2,6,100);
set(gca,'YScale','linear','XScale','linear','XLim',xAxisBootRange,'YLim',yAxisRange);
xlabel(xLabelBootStr,'Interpreter','latex'); 
ylabel(yLabelBootStr,'Interpreter','latex'); 
text(0.6,0.5,['\sigma_{single-onset} = ',num2str(f1.p1,4),' \Omega^{-1} m^{-1}'],'Units','Normalized','Color',cmap(60,:));
text(0.6,0.4,['\sigma_{multi-onset} = ',num2str(f2.p1,4),' \Omega^{-1} m^{-1}'],'Units','Normalized','Color',cmap(200,:));
legend([p1,p2],{'Single-onset','Multi-onset'});

%% CCA of substorm types with 

%% Plot of pdf(epsilon,sml)

clearvars multiStorm singleStorm neitherStorm
SML = -100:-100:-900;

for j=1:length(SML)-1
    SML_max = SML(j);
    SML_min = SML(j+1);
    multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0) & T.peakSML<SML_max & T.peakSML>SML_min;
    singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0) & T.peakSML<SML_max & T.peakSML>SML_min;
    % multiSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration<duration(3,0,0);
    % singleSubstorm = T.previousSubstormDuration>duration(3,0,0) & T.nextSubstormDuration>duration(3,0,0);
    neitherSubstorm = ~(multiSubstorm | singleSubstorm) & T.peakSML<SML_max & T.peakSML>SML_min;

    multiStormID = T.stormID(multiSubstorm);
    singleStormID = T.stormID(singleSubstorm);
    neitherStormID = T.stormID(neitherSubstorm);
    sampleSize = 30;
    ensembleSize = 5000;
    for i= 1:ensembleSize
        multiStorm.EnsembleID(j,i,:)=multiStormID(randperm(length(multiStormID),sampleSize));
        singleStorm.EnsembleID(j,i,:)=singleStormID(randperm(length(singleStormID),sampleSize));
        neitherStorm.EnsembleID(j,i,:)=neitherStormID(randperm(length(neitherStormID),sampleSize));
    end
end

%%
EArray1 = T.EArray; 
SMLArray1 = T.smlArray;
MLATArray1 = T.MLAT;

multiStorm.EArray = squeeze(nanmean(reshape(EArray1(multiStorm.EnsembleID,:),size(multiStorm.EnsembleID,1),size(multiStorm.EnsembleID,2),size(multiStorm.EnsembleID,3),size(EArray1,2)),3));
singleStorm.EArray = squeeze(nanmean(reshape(EArray1(singleStorm.EnsembleID,:),size(singleStorm.EnsembleID,1),size(singleStorm.EnsembleID,2),size(singleStorm.EnsembleID,3),size(EArray1,2)),3));
neitherStorm.EArray = squeeze(nanmean(reshape(EArray1(neitherStorm.EnsembleID,:),size(neitherStorm.EnsembleID,1),size(neitherStorm.EnsembleID,2),size(neitherStorm.EnsembleID,3),size(EArray1,2)),3));

multiStorm.SMLArray = squeeze(nanmean(reshape(SMLArray1(multiStorm.EnsembleID,:),size(multiStorm.EnsembleID,1),size(multiStorm.EnsembleID,2),size(multiStorm.EnsembleID,3),size(SMLArray1,2)),3));
singleStorm.SMLArray = squeeze(nanmean(reshape(SMLArray1(singleStorm.EnsembleID,:),size(singleStorm.EnsembleID,1),size(singleStorm.EnsembleID,2),size(singleStorm.EnsembleID,3),size(SMLArray1,2)),3));
neitherStorm.SMLArray = squeeze(nanmean(reshape(SMLArray1(neitherStorm.EnsembleID,:),size(neitherStorm.EnsembleID,1),size(neitherStorm.EnsembleID,2),size(neitherStorm.EnsembleID,3),size(SMLArray1,2)),3));

multiStorm.MLATArray = squeeze(nanmean(reshape(MLATArray1(multiStorm.EnsembleID,:),size(multiStorm.EnsembleID,1),size(multiStorm.EnsembleID,2),size(multiStorm.EnsembleID,3),size(MLATArray1,2)),3));
singleStorm.MLATArray = squeeze(nanmean(reshape(MLATArray1(singleStorm.EnsembleID,:),size(singleStorm.EnsembleID,1),size(singleStorm.EnsembleID,2),size(singleStorm.EnsembleID,3),size(MLATArray1,2)),3));
neitherStorm.MLATArray = squeeze(nanmean(reshape(MLATArray1(neitherStorm.EnsembleID,:),size(neitherStorm.EnsembleID,1),size(neitherStorm.EnsembleID,2),size(neitherStorm.EnsembleID,3),size(MLATArray1,2)),3));



singleStorm.timeArrayRelOnset = T.timeArrayRelOnset(1,:);

%%
% growthTIndex = 60:180; %Minutes
% expansionTIndex = 181:301; %Minutes

growthTIndex = 1:2; %Minutes
expansionTIndex = 360:361; %Minutes
onsetTIndex = 160:200;

% plot histograms of Earray for specific AE

figure; 
subplot(2,1,1)
% edges = linspace(0.5e+4,4.5e+4,50);
edges = linspace(0.5e+1,4.5e+2,50);
histogram(trapz(singleStorm.EArray(2,:,growthTIndex),3),edges,'Normalization','countdensity'); hold on; 
histogram(trapz(multiStorm.EArray(2,:,growthTIndex),3),edges,'Normalization','countdensity');  
% histogram(mean(singleStorm.EArray(2,:,growthTIndex),3),edges,'Normalization','countdensity'); hold on; 
% histogram(mean(multiStorm.EArray(2,:,growthTIndex),3),edges,'Normalization','countdensity'); 
set(gca,'XScale','linear');
text(0.35,0.8,[num2str(SML(2)),'<SML<',num2str(SML(3)),'nT'],'Units','Normalized');
xlabel('$\int_{t_0-2\mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
legend(['Single-onset'],...
    ['Multi-onset' ]);
% xlim([0.5e+4,4.5e+4]);
set(gca,'XScale','linear');
title('Pre-Onset');
ylabel('Count Density');

subplot(2,1,2)
histogram(trapz(singleStorm.EArray(2,:,expansionTIndex),3),edges,'Normalization','countdensity'); hold on; 
histogram(trapz(multiStorm.EArray(2,:,expansionTIndex),3),edges,'Normalization','countdensity');  
% histogram(mean(singleStorm.EArray(2,:,expansionTIndex),3),edges,'Normalization','countdensity'); hold on; 
% histogram(mean(multiStorm.EArray(2,:,expansionTIndex),3),edges,'Normalization','countdensity');  
set(gca,'XScale','linear');
text(0.35,0.8,[num2str(SML(2)),'<SML<',num2str(SML(3)),'nT'],'Units','Normalized');
xlabel('$\int_{t_0}^{t_0+2\mathrm{Hr}} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
legend(['Single-onset'],...
    ['Multi-onset']);
% xlim([0.5e+4,4.5e+4]);
set(gca,'XScale','linear');
title('Post-Onset');
ylabel('Count Density');
%% Plot SML, E distribution with contour plots

fig=figure; 
resize_figure(fig,200,250);
edges = linspace(2e+5,70e+5,200);
smledges = fliplr(-100:-100:-1000);
% subplot(2,1,1)


singleStorm.EShaped = reshape(singleStorm.EArray,[size(singleStorm.EArray,1)*size(singleStorm.EArray,2),size(singleStorm.EArray,3)]);
singleStorm.SMLShaped = reshape(singleStorm.SMLArray,[size(singleStorm.SMLArray,1)*size(singleStorm.SMLArray,2),size(singleStorm.SMLArray,3)]);
singleStorm.MLATShaped = reshape(singleStorm.MLATArray,[size(singleStorm.MLATArray,1)*size(singleStorm.MLATArray,2),1]);

multiStorm.EShaped = reshape(multiStorm.EArray,[size(multiStorm.EArray,1)*size(multiStorm.EArray,2),size(multiStorm.EArray,3)]);
multiStorm.SMLShaped = reshape(multiStorm.SMLArray,[size(multiStorm.SMLArray,1)*size(multiStorm.SMLArray,2),size(multiStorm.SMLArray,3)]);
multiStorm.MLATShaped = reshape(multiStorm.MLATArray,[size(multiStorm.MLATArray,1)*size(multiStorm.MLATArray,2),1]);

neitherStorm.EShaped = reshape(neitherStorm.EArray,[size(neitherStorm.EArray,1)*size(neitherStorm.EArray,2),size(neitherStorm.EArray,3)]);
neitherStorm.SMLShaped = reshape(neitherStorm.SMLArray,[size(neitherStorm.SMLArray,1)*size(neitherStorm.SMLArray,2),size(neitherStorm.SMLArray,3)]);
multiStorm.MLATShaped = reshape(multiStorm.MLATArray,[size(multiStorm.MLATArray,1)*size(multiStorm.MLATArray,2),1]);

[h1.Value, h1.X, h1.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),...
    singleStorm.EShaped(:,growthTIndex),2),min(singleStorm.SMLShaped(:,onsetTIndex),[],2),...
    edges, smledges,...
    'Normalization','countdensity'); 

    [H1X, H1Y]=meshgrid(0.5.*(h1.X(1:end-1)+h1.X(2:end)),(h1.Y(1:end-1)+h1.Y(2:end))*0.5);
    
    axes1 = axes;


    
    [C,Ch]=contourf(H1X,H1Y,h1.Value',logspace(-4,-6,4));
    colormap(axes1,get_colormap([1,1,1],[0,0,0]));
    axes1.XScale = 'linear';
    
    set(axes1,'ColorScale','log','XLim',[0.1e+6, 8e+6]);
    
    caxis(axes1,[1e-6,1e-4]);
    
    c1 = colorbar_thin('Location','eastoutside','YLabel',{'Single-Onset Substorms','[Count Density]'});
%     clabel(C,Ch,[1e-6,1e-4],'Color',[0.5,0.5,0.5]);
    xlabel('$\int_{t_0-2 \mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
    title('Pre-Onset');
    ylabel('Peak SML [nT]');
    zlabel('Count Density');
    
    hold on;
    
[h2.Value, h2.X, h2.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),...
    multiStorm.EShaped(:,growthTIndex),2),min(multiStorm.SMLShaped(:,onsetTIndex),[],2),...
    edges, smledges,...
    'Normalization','countdensity');
    
    [H2X, H2Y]=meshgrid(0.5.*(h2.X(1:end-1)+h2.X(2:end)),(h2.Y(1:end-1)+h2.Y(2:end))*0.5);
    
    axes2 = axes;
    
    
    [C2,Ch2]=contour(H2X,H2Y,h2.Value',logspace(-4,-6,4));
    axes2.Visible= 'off';
    axes2.XScale = 'linear';
    alpha 0.9;
    c2 = colorbar_thin('Location','westoutside','YLabel',{'Multi-Onset Substorms [Count Density]'});
%     clabel(C2,Ch2,[1e-6,1e-4],'Color','r');
    set(axes2,'ColorScale','log','XLim',[0.1e+6, 8e+6]);
    colormap(axes2,get_colormap('m','r'));
    caxis(axes2,[1e-6,1e-4]);
    linkaxes([axes1, axes2]);

    text(0.01,0.3,{['#Substorms-per-Ensemble: ',num2str(sampleSize)],['#Ensembles: ',num2str(ensembleSize)]},'Units','Normalized');
    
    
fig2=figure; 
resize_figure(fig2,200,250);
[h1.Value, h1.X, h1.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),...
    singleStorm.EShaped(:,expansionTIndex),2),min(singleStorm.SMLShaped(:,onsetTIndex),[],2),...
    edges, smledges,...
    'Normalization','countdensity'); 

    [H1X, H1Y]=meshgrid(0.5.*(h1.X(1:end-1)+h1.X(2:end)),(h1.Y(1:end-1)+h1.Y(2:end))*0.5);
    
    axes1 = axes;


    
    contourf(H1X,H1Y,h1.Value',logspace(-4,-6,4));
    colormap(axes1,get_colormap([1,1,1],[0,0,0]));
    axes1.XScale = 'linear';
    
    set(axes1,'ColorScale','log','XLim',[0.1e+6, 8e+6]);
    
    caxis(axes1,[1e-6,1e-4]);
    
    c1 = colorbar_thin('Location','eastoutside','YLabel',{'Single-Onset Substorms','[Count Density]'});
    
    xlabel('$\int_{t_0 \mathrm{Hr}}^{t_0+2} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
    title('Post-Onset');
    ylabel('Peak SML [nT]');
    
    hold on;
    
[h2.Value, h2.X, h2.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),...
    multiStorm.EShaped(:,expansionTIndex),2),min(multiStorm.SMLShaped(:,onsetTIndex),[],2),...
    edges, smledges,...
    'Normalization','countdensity');
    
    [H2X, H2Y]=meshgrid(0.5.*(h2.X(1:end-1)+h2.X(2:end)),(h2.Y(1:end-1)+h2.Y(2:end))*0.5);
    
    axes2 = axes;
    
    
    contour(H2X,H2Y,h2.Value',logspace(-4,-6,4));
    axes2.Visible= 'off';
    axes2.XScale = 'linear';
    alpha 0.9;
    c2 = colorbar_thin('Location','westoutside','YLabel',{'Multi-Onset Substorms [Count Density]'});
    set(axes2,'ColorScale','log','XLim',[0.1e+6, 8e+6]);
    colormap(axes2,get_colormap('m','r'));
    caxis(axes2,[1e-6,1e-4]);
    linkaxes([axes1, axes2]);

    text(0.01,0.3,{['#Substorms-per-Ensemble: ',num2str(sampleSize)],['#Ensembles: ',num2str(ensembleSize)]},'Units','Normalized');

    %% Plot SML, E, MLAT distribution with contour plots

fig=figure; 
resize_figure(fig,200,250);
edges = logspace(5,7,200);
smledges = fliplr(-100:-100:-1000);
mlatedges = 60:1:80;
% subplot(2,1,1)


singleStorm.EShaped = reshape(singleStorm.EArray,[size(singleStorm.EArray,1)*size(singleStorm.EArray,2),size(singleStorm.EArray,3)]);
singleStorm.SMLShaped = reshape(singleStorm.SMLArray,[size(singleStorm.SMLArray,1)*size(singleStorm.SMLArray,2),size(singleStorm.SMLArray,3)]);
singleStorm.MLATShaped = reshape(singleStorm.MLATArray,[size(singleStorm.MLATArray,1)*size(singleStorm.MLATArray,2),1]);

multiStorm.EShaped = reshape(multiStorm.EArray,[size(multiStorm.EArray,1)*size(multiStorm.EArray,2),size(multiStorm.EArray,3)]);
multiStorm.SMLShaped = reshape(multiStorm.SMLArray,[size(multiStorm.SMLArray,1)*size(multiStorm.SMLArray,2),size(multiStorm.SMLArray,3)]);
multiStorm.MLATShaped = reshape(multiStorm.MLATArray,[size(multiStorm.MLATArray,1)*size(multiStorm.MLATArray,2),1]);

neitherStorm.EShaped = reshape(neitherStorm.EArray,[size(neitherStorm.EArray,1)*size(neitherStorm.EArray,2),size(neitherStorm.EArray,3)]);
neitherStorm.SMLShaped = reshape(neitherStorm.SMLArray,[size(neitherStorm.SMLArray,1)*size(neitherStorm.SMLArray,2),size(neitherStorm.SMLArray,3)]);
multiStorm.MLATShaped = reshape(multiStorm.MLATArray,[size(multiStorm.MLATArray,1)*size(multiStorm.MLATArray,2),1]);

[h1.Value, h1.X, h1.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),...
    singleStorm.EShaped(:,growthTIndex),2),singleStorm.MLATShaped,...
    edges, mlatedges,...
    'Normalization','countdensity'); 

    [H1X, H1Y]=meshgrid(0.5.*(h1.X(1:end-1)+h1.X(2:end)),(h1.Y(1:end-1)+h1.Y(2:end))*0.5);
    
    axes1 = axes;


    
    [C,Ch]=contourf(H1X,H1Y,h1.Value',logspace(-2,-8,5));
    colormap(axes1,get_colormap([1,1,1],[0,0,0]));
    axes1.XScale = 'log';
    
    set(axes1,'ColorScale','log','XLim',[10^5, 10^7]);
    
    caxis(axes1,[1e-6,1e-2]);
    
    c1 = colorbar_thin('Location','eastoutside','YLabel',{'Single-Onset Substorms','[Count Density]'});
%     clabel(C,Ch,[1e-6,1e-4],'Color',[0.5,0.5,0.5]);
    xlabel('$\int_{t_0-2 \mathrm{Hr}}^{t_0} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
    title('Pre-Onset');
    ylabel('MLAT [^\circ N]');
    zlabel('Count Density');
    
    hold on;
    
[h2.Value, h2.X, h2.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,growthTIndex)),...
    multiStorm.EShaped(:,growthTIndex),2),multiStorm.MLATShaped,...
    edges, mlatedges,...
    'Normalization','countdensity');
    
    [H2X, H2Y]=meshgrid(0.5.*(h2.X(1:end-1)+h2.X(2:end)),(h2.Y(1:end-1)+h2.Y(2:end))*0.5);
    
    axes2 = axes;
    
    
    [C2,Ch2]=contour(H2X,H2Y,h2.Value',logspace(-2,-8,5));
    axes2.Visible= 'off';
    axes2.XScale = 'log';
    alpha 0.9;
    c2 = colorbar_thin('Location','westoutside','YLabel',{'Multi-Onset Substorms [Count Density]'});
%     clabel(C2,Ch2,[1e-6,1e-4],'Color','r');
    set(axes2,'ColorScale','log','XLim',[10^5, 10^7]);
    colormap(axes2,get_colormap('m','r'));
%     caxis(axes2,[1e-6,1e-4]);
    linkaxes([axes1, axes2]);

    text(0.01,0.3,{['#Substorms-per-Ensemble: ',num2str(sampleSize)],['#Ensembles: ',num2str(ensembleSize)]},'Units','Normalized');
    

fig2=figure; 
resize_figure(fig2,200,250);
[h1.Value, h1.X, h1.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),...
    singleStorm.EShaped(:,expansionTIndex),2),singleStorm.MLATShaped,...
    edges, mlatedges,...
    'Normalization','countdensity'); 

    [H1X, H1Y]=meshgrid(0.5.*(h1.X(1:end-1)+h1.X(2:end)),(h1.Y(1:end-1)+h1.Y(2:end))*0.5);
    
    axes1 = axes;


    
    contourf(H1X,H1Y,h1.Value',logspace(-2,-8,5));
    colormap(axes1,get_colormap([1,1,1],[0,0,0]));
    axes1.XScale = 'log';
    
    set(axes1,'ColorScale','log','XLim',[10^5, 10^7]);
    
    caxis(axes1,[1e-6,1e-2]);
    
    c1 = colorbar_thin('Location','eastoutside','YLabel',{'Single-Onset Substorms','[Count Density]'});
    
    xlabel('$\int_{t_0 \mathrm{Hr}}^{t_0+2} V B^2 \sin^4 \frac{\theta_c}{2} L_0^2 \ \mathrm{d}t$ [GJ]','Interpreter','latex'); 
    title('Post-Onset');
    ylabel('MLAT [^\circ N]');
    
    hold on;
    
[h2.Value, h2.X, h2.Y] = histcounts2(trapz(seconds(T.timeArrayRelOnset(1,expansionTIndex)),...
    multiStorm.EShaped(:,expansionTIndex),2),multiStorm.MLATShaped,...
    edges, mlatedges,...
    'Normalization','countdensity');
    
    [H2X, H2Y]=meshgrid(0.5.*(h2.X(1:end-1)+h2.X(2:end)),(h2.Y(1:end-1)+h2.Y(2:end))*0.5);
    
    axes2 = axes;
    
    
    contour(H2X,H2Y,h2.Value',logspace(-2,-8,5));
    axes2.Visible= 'off';
    axes2.XScale = 'log';
    alpha 0.9;
    c2 = colorbar_thin('Location','westoutside','YLabel',{'Multi-Onset Substorms [Count Density]'});
    set(axes2,'ColorScale','log','XLim',[10^5, 10^7]);
    colormap(axes2,get_colormap('m','r'));
    caxis(axes2,[1e-6,1e-2]);
    linkaxes([axes1, axes2]);

    text(0.01,0.3,{['#Substorms-per-Ensemble: ',num2str(sampleSize)],['#Ensembles: ',num2str(ensembleSize)]},'Units','Normalized');


%%
plot_superposed_epoch(T, logical(sum(T.stormID==singleStorm.EnsembleID(1964,:),2)) , {'BzArray','EArray','symHArray'});
%%
plot_epoch(T, multiStorm.EnsembleID(1964,3 ), {'BzArray','EArray'});


%% Functions

prevSuperMagFileType = superMagFileType; 


function plot_superposed_epoch(T, condition, parameters)
    
    tRange = [min(T.timeArrayRelOnset(1,:)) max(T.timeArrayRelOnset(1,:))];

    h=figure; 
    totalNo = length(parameters);
    p=create_panels(h,'totalPanelNo',totalNo,'marginbottom',10,'panelHeight',40);
    
    for i = 1:1:totalNo
        p(1,i).select();
        plot_specified_parameter(T,condition,tRange,parameters{i});
        
        if i~=totalNo
            set(gca,'XTick',{});    
        else
            xlabel('Substorm Time (t)');
        end
    end

end

function plot_epoch(T, stormID, parameters)
    
    tRange = [min(T.timeArrayRelOnset(1,:)) max(T.timeArrayRelOnset(1,:))];

    h=figure; 
    totalNo = length(parameters);
    p=create_panels(h,'totalPanelNo',totalNo,'marginbottom',10,'panelHeight',40);
    
    for i = 1:1:totalNo
        p(1,i).select();
        plot_specified_parameter_one_storm(T,stormID,tRange,parameters{i});
        
        if i~=totalNo
            set(gca,'XTick',{});    
        else
            xlabel('Substorm Time (t)');
        end
    end

end

function plot_specified_parameter(T,condition,tRange,parameterString)
    
    switch(parameterString)
        case 'BzArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'BzArray'});
            ylabel('IMF Bz [nT]')
            yRange = [-3 1];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'BzArray'}));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],yRange,[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'thetaArray'
            time_plot(T.timeArrayRelOnset(1,:),sin(T{condition,'thetaArray'}./2).^2);
            ylabel('sin^2(\Theta/2)');
            yRange = [0.4 0.8];
            ylim(yRange);
            med = nanmedian(nanmedian(sin(T{condition,'thetaArray'}./2).^2));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],[0 1],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'pArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'pArray'});
            ylabel('P_{dyn} [nPa]');
            yRange = [1.6,2.8];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'pArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[1.6 2.8],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'smlArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'smlArray'});
            ylabel('SML [nT]')
            yRange = [-400,0];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'smlArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-400 0],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
         
        case 'symHArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'symHArray'},15);
            ylabel('SYM-H [nT]')
            yRange = [-30,10];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'symHArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-30 10],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'r_0Array' 
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'r_0Array'});
            ylabel('r_0 [R_E]')
            yRange = [9.5,10.5];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'r_0Array'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[9 11],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'alphaArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'alphaArray'});
            ylabel('\alpha [a.u.]')
            % ylim([9.5,10.5]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');

        case 'EArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'EArray'});
            ylabel('\epsilon(t) [GW]');
            hold on;
            med = nanmedian(nanmedian(T{condition,'BzArray'}));
            plot3(tRange,[med med] ,[1 1]);
            hold on;
            plot3([duration duration],[0 600],[1 1]);
            % set(gca,'YScale','log');
            ylim([10 400]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');           
            
        case 'pfisrNe'
            text(0.5,0.5,1,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation))))],'Units','normalized');
           
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrNe'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^10 10^12]);
            ylabel('N_e');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        case 'pfisrEflux'
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrEflux'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrEbins'}{1});
            set(gca,'ColorScale','log','CLim',[10^8 10^12],'YScale','log');
            xlim([-3,3]);
            ylabel('\phi(E)');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        case 'pfisrHall'
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrHall'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^-6 10^-2]);
            ylabel('\Sigma_H');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        otherwise
            error(['Invalid parameter ',parameterString]);
    end
    
end

function plot_specified_parameter_one_storm(T,stormID,tRange,parameterString)
    
    switch(parameterString)
        case 'BzArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'BzArray'});
            ylabel('IMF Bz [nT]')
            yRange = [-3 1];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'BzArray'}));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],yRange,[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'thetaArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),sin(T{stormID,'thetaArray'}./2).^2);
            ylabel('sin^2(\Theta/2)');
            yRange = [0.4 0.8];
            ylim(yRange);
            med = nanmedian(nanmedian(sin(T{stormID,'thetaArray'}./2).^2));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],[0 1],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'pArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'pArray'});
            ylabel('P_{dyn} [nPa]');
            yRange = [1.6,2.8];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'pArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[1.6 2.8],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'smlArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'smlArray'});
            ylabel('SML [nT]')
            yRange = [-400,0];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'smlArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-400 0],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
         
        case 'symHArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'symHArray'});
            ylabel('SYM-H [nT]')
            yRange = [-30,10];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'symHArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-30 10],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'r_0Array' 
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'r_0Array'});
            ylabel('r_0 [R_E]')
            yRange = [9.5,10.5];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'r_0Array'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[9 11],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'alphaArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'alphaArray'});
            ylabel('\alpha [a.u.]')
            % ylim([9.5,10.5]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');

        case 'EArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'EArray'});
            ylabel('\epsilon(t) [GW]');
            hold on;
            med = nanmedian(nanmedian(T{stormID,'BzArray'}));
            plot3(tRange,[med med] ,[1 1]);
            hold on;
            plot3([duration duration],[0 600],[1 1]);
            % set(gca,'YScale','log');
            ylim([10 400]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'pfisrNe'
            if ~ismissing(T.storageLocation(stormID))
                time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                    T{stormID ,'pfisrNe'}{1},...
                    T{stormID ,'pfisrAlt'}{1});
                set(gca,'ColorScale','log','CLim',[10^10 10^12]);
                ylabel('N_e');
                text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
                text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
            
        case 'pfisrEflux'
            if ~ismissing(T.storageLocation(stormID))
            time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                T{stormID ,'pfisrEflux'}{1},...
                T{stormID ,'pfisrEbins'}{1});
            set(gca,'ColorScale','log','CLim',[10^8 10^12],'YScale','log');
            xlim([-3,3]);
            ylabel('\phi(E)');
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
            
        case 'pfisrHall'
            if ~ismissing(T.storageLocation(stormID))
            time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                T{stormID ,'pfisrHall'}{1},...
                T{stormID ,'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^-6 10^-2]);
            ylabel('\Sigma_H');
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
        otherwise
            error(['Invalid parameter ',parameterString]);
    end
    
end

function time_plot_stormID(t,y)

plot3(t,y,20.*ones(1,numel(t)),'r');
view(2);

end

function time_plot(t,y, pdfBinNum)

if nargin<3
    pdfBinNum = 300;
end

% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;
my = nanmean(y);
medy = nanmedian(y);
sdev = nanstd(y);
sy = sdev./((size(y,1).^0.5));

edges = linspace(min(my)-nanmean(sdev),max(my)+nanmean(sdev),pdfBinNum);
for i=1:size(y,2) 
    N(:,i) = histcounts(y(:,i),edges,'Normalization','pdf');
end
[T,Y] = meshgrid(t,edges(1:end-1));


s=surf(T,Y,N); 
view(2); 
s.EdgeColor = 'none'; 
colormap(get_colormap('w','k'));
colorbar_thin('YLabel','PDF');
hold on;


plot3(t,my,20.*ones(1,numel(t)),'r'); hold on;
plot3(t,medy,20.*ones(1,numel(t)),'m'); hold on;
plot3(t,my+sy,20.*ones(1,numel(t)),'k'); hold on;
plot3(t,my-sy,20.*ones(1,numel(t)),'k');



end

function time_plot_stormID_2D(t,y,z)


% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;

[X,Y] = meshgrid(t,z);

zValue = y';
    
h=pcolor(hours(X),Y,zValue);
set(h,'EdgeColor','none');
shading flat;
set(gca,'Layer','top');

hold on;

%% Calculating Y axis and Color axis limits
    
	miny=min(z); maxy=max(z);
	ylim([miny maxy]);
    

    
    minz=min(min(zValue)); maxz=max(max(zValue));
    if ~isnan(minz) && ~isnan(maxz)
        caxis([minz maxz]);
    end
    colormap(get_colormap('w','k'));
    colorbar_thin();

end


function time_plot_2D(t,y,z)


% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;
yarr = cat(3,y{:});
my = nanmean(yarr,3);
medy = nanmedian(yarr,3);
sdev = squeeze(nanstd(permute(yarr,[3 1 2])));
sy = sdev./((size(yarr,1).^0.5));

[X,Y] = meshgrid(t,z);

zValue = my';
    
h=pcolor(hours(X),Y,zValue);
set(h,'EdgeColor','none');
shading flat;
set(gca,'Layer','top');

hold on;

%% Calculating Y axis and Color axis limits
    
	miny=min(z); maxy=max(z);
	ylim([miny maxy]);
    

    
    minz=min(min(zValue)); maxz=max(max(zValue));
    if ~isnan(minz) && ~isnan(maxz)
        caxis([minz maxz]);
    end
    colormap(get_colormap('w','k'));
    colorbar_thin();

end


function T = add_pfisr_array_to_table(T,alt)
    % Adds pfisr data, electron density, Hall, and energy flux on to the
    % table
    
    if nargin<2
        alt = 60:1:140;
    end
    
    for i = T.stormID(~ismissing(T.storageLocation))'
        try
        pfisrNe = h5read(T{i,'storageLocation'},'/inputData/Ne');
        pfisrHall = h5read(T{i,'storageLocation'},'/conductivity/hall');
        pfisrEflux = h5read(T{i,'storageLocation'},'/energy/energyFlux');
        pfisrEbins = h5read(T{i,'storageLocation'},'/energy/energyBin'); 
        pfisrtime = h5read(T{i,'storageLocation'},'/time');
        pfisralt =  h5read(T{i,'storageLocation'},'/alt');
        pfisrRelTime = interp1(T.timeArray(i,:),T.timeArrayRelOnset(i,:),pfisrtime,'linear','extrap');
        
        T.pfisrNe{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrNe,alt,'linear')',T.timeArrayRelOnset(i,:));
        T.pfisrHall{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrHall,alt,'linear')',T.timeArrayRelOnset(i,:));
        T.pfisrEbins{i} = pfisrEbins;
        T.pfisrEflux{i} = interp1(pfisrRelTime,pfisrEflux',T.timeArrayRelOnset(i,:));
        T.pfisrAlt{i} = alt;
        catch ME
             disp(['Substorm ID: ',num2str(i)]);
             getReport(ME)
             warning('Continuing despite error');
        end
    end

end


function T = add_omni_array_to_table(T, omni, preOnsetWindow, postOnsetWindow)

    if nargin<4
        postOnsetWindow = duration(3,0,0);
    end
    
    if nargin<3
        preOnsetWindow = duration(3,0,0);
    end
    
    smlTime = -preOnsetWindow:duration(0,1,0):postOnsetWindow;
    T.timeArray = repmat(duration,height(T),length(smlTime));
    zeroMatrix = zeros(height(T),length(smlTime));
    T.smlArray = zeroMatrix;
    T.smuArray = zeroMatrix;
    T.ALArray = zeroMatrix;
    T.BzArray = zeroMatrix;
    T.pArray = zeroMatrix;
    T.densityArray = zeroMatrix;
    T.vArray = zeroMatrix;
    T.EklArray = zeroMatrix;
    T.EArray = zeroMatrix;
    T.thetaArray = zeroMatrix;
    T.r_0Array = zeroMatrix;
    T.alphaArray = zeroMatrix;
    T.symHArray = zeroMatrix;
    
    for i=1:1:height(T)
    
        TMatrix(i,:) = datenum(T.Time(i)-preOnsetWindow : duration(0,1,0) : ...
        T.Time(i)+postOnsetWindow);
    
    end
    
    T.timeArrayRelOnset = repmat(smlTime,height(T),1);
    T.timeArray = TMatrix;
    
    T.smlArray      = omni.Fsml(TMatrix);
    T.smuArray      = omni.Fsmu(TMatrix);
    T.pArray        = omni.Fp(TMatrix); 
    T.BzArray       = omni.FBz(TMatrix); 
    T.ALArray       = omni.FAL(TMatrix); 
    T.densityArray  = omni.Fdensity(TMatrix); 
    T.vArray        = omni.Fv(TMatrix); 
    T.eklArray      = omni.Fekl(TMatrix); 
    T.EArray        = omni.FE(TMatrix); 
    T.thetaArray    = omni.Ftheta(TMatrix); 
    T.r_0Array      = omni.Fr_0(TMatrix); 
    T.alphaArray    = omni.Falpha(TMatrix);
    T.symHArray     = omni.FsymH(TMatrix);

end



function T = identify_qualified_substorm_duration(T)
    % Function to identify and define isolated substorms
    % identify_qualified_substorm_duration, 
    % first filters substorms from the database that are likely legitimate
    % removing activity that occurs in the day-side
    % secondly calculates the duration between the above qualified
    % substorms
    % Adds two columns to Table T
    %   previousSubstormDuration : The amount of time before which the
    %                              previous substorm had its onset. 
    %   nextSubstormDuration     : The amount of time after which the
    %                              next substorm will have its onset
    
    MLT_qualification = T.MLT >=16 | T.MLT<=8; 
    SML_qualification = T.peakSML <= -100; 
    
    qualified_substorms = MLT_qualification & SML_qualification; 
    T.qualifiedSubstorms = qualified_substorms;
    
    stormID = T.stormID(qualified_substorms);
    time = T.Time(qualified_substorms);
    substormDuration = repmat(duration, length(time),1);
    substormDuration(2:end) = time(2:end) - time(1:end-1);
    substormDuration(1) = nan;
    
    nextSubstormDuration = repmat(duration, length(time),1);
    nextSubstormDuration(1:end-1) = time(2:end) - time(1:end-1);
    nextSubstormDuration(end) = nan;
    
    
    T.previousSubstormDuration = repmat(duration, length(T.Time),1);
    T.previousSubstormDuration(:) = nan(length(T.Time),1);
    T.previousSubstormDuration(T.stormID(stormID)) = substormDuration;
    
    T.nextSubstormDuration = repmat(duration, length(T.Time),1);
    T.nextSubstormDuration(:) = nan(length(T.Time),1);
    T.nextSubstormDuration(T.stormID(stormID)) = nextSubstormDuration;
    
end

% Extracting Omni data and developing gridded interpolants of the database
function omni = extract_omni_data(omniFile)
    
    %Time
    time = unixtime2matlab(h5read(omniFile,'/Time'));
    
        % Auroral electrojet indices
    SML = h5read(omniFile,'/Indices/SML');
    omni.Fsml = griddedInterpolant(time,SML);
    
    SMU = h5read(omniFile,'/Indices/SMU');
    omni.Fsmu = griddedInterpolant(time,SMU);
    
    AL = h5read(omniFile,'/Indices/AL');
    AL(AL==99999)=nan;
    omni.FAL = griddedInterpolant(time, AL);
    
    symH = h5read(omniFile,'/Indices/SYM_H');
    omni.FsymH = griddedInterpolant(time, symH);
    % Solar wind
    
        % Dynamic Pressure
    Pdyn = h5read(omniFile,'/FlowPressure');
    omni.Fp = griddedInterpolant(time,Pdyn);
    
        % IMF Bz
    BzGSM = h5read(omniFile,'/BField/BzGSM');
    omni.FBz = griddedInterpolant(time,BzGSM);
    
        % IMF By
    ByGSM = h5read(omniFile,'/BField/ByGSM');
    omni.FBy = griddedInterpolant(time,ByGSM);
    
        % IMF B magnitude
    BxGSE = h5read(omniFile,'/BField/BxGSE');
    ByGSE = h5read(omniFile,'/BField/ByGSE');
    BzGSE = h5read(omniFile,'/BField/BzGSE');
    B = (BxGSE.^2+ByGSE.^2+BzGSE.^2).^0.5;
    omni.FB = griddedInterpolant(time,B);
    
        % IMF B_T (Tangential to to GSM_x)
    BzGSM(BzGSM==0) = 0.0001;
    B_T = (ByGSM.^2 + BzGSM.^2).^0.5;
    
        % IMF Clock angle
    theta_c = wrapTo2Pi(atan2(ByGSM,BzGSM));
    theta_kl = wrapTo2Pi(atan2(B_T,BzGSM));
    omni.Ftheta = griddedInterpolant(time,theta_c);
        
        % Density
    density = h5read(omniFile,'/ProtonDensity');
    omni.Fdensity = griddedInterpolant(time,density);
        
        % Velocity
    velocity = h5read(omniFile,'/Velocity/V');
    omni.Fv = griddedInterpolant(time,velocity);
    
    
    % Solarwind - Magnetosphere Coupling
    l_0 = 7*(6371*10^3);
    
    E_kl = velocity.*B_T.*(sin(theta_kl/2)).^2; %Km nT/s
    % The “geoeffective” (or “merging”) electric field [Kan and Lee, 1979]
    omni.Fekl = griddedInterpolant(time,E_kl);
    
    E = 1e-9.*1e7.*(velocity.*10.^3).*((B.*10^-9).^2).*(sin(theta_c/2)).^4*l_0^2; %GW 
    omni.FE = griddedInterpolant(time,E);
    
    
    % Magnetopause parameters from Shen et al., 1993

        % r_0 is the standoff distance in R_E
    r_0 = (11.4 + 0.013.*BzGSM).*(Pdyn).^(-1./6.6); % for Bz>=0
    r_0(BzGSM<0) = (11.4 + 0.14.*BzGSM(BzGSM<0)).*(Pdyn(BzGSM<0)).^(-1./6.6); % for Bz<0
    omni.Fr_0 = griddedInterpolant(time,r_0);
    
        % alpha is the level of tail flaring
    alpha = (0.58-0.010*BzGSM).*(1+0.01*Pdyn);
    omni.Falpha = griddedInterpolant(time,alpha);
    
    
end


% Creating a function that does what substorm_table.m does, i.e., outputs a
% table of substorms that are linked with PFISR and DASC data given an
% input substorm file. 

function [T, superMag] = substorm_create_table(superMagFileStr, superMagFileType,...
  timeMinStr, timeMaxStr)

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/Elements/Nithin/Data/';
    storeDir = '/media/nithin/Elements/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
else
%     error(['Not configured for this computer: ',get_computer_name]);
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
end
outputAMISRFileStr = 'amisrWebDatabase.h5';
amisrDatabaseStr = [dataDir,outputAMISRFileStr];
% dascFileStr = [storeDir,'dascDatabase.h5'];
omniFileStr = [dataDir,'omni.h5'];

% Substorms at PFISR [IMPORTANT] : Range of closeness to PFISR
% Dmlt = 2;
% Dmlat = 50; %desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

% Loading database

% Create amisr database
if ~isfile(amisrDatabaseStr)
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataDir,outputAMISRFileStr]);
end
%
% Load supermag database
superMag = extract_superMag_data(superMagFileStr,superMagFileType);
superMag.stormID = (1:length(superMag.time))';

% Load amisr data
amisr = extract_amisr_data(amisrDatabaseStr);

% Load omni data
omni.AE = h5read(omniFileStr,'/Indices/AE');
omni.time = unixtime2matlab(h5read(omniFileStr,'/Time'));

% Calculations
% Estimating PFISR magnetic coordinates 
[superMag.pfisrMlat,superMag.pfisrMlon,superMag.pfisrMlt] = get_magnetic_coordinates([pkrh0,pkrGLAT,pkrGLON],superMag.time(:));

% Interpolating AE index of the substorms
superMag.AE = interp1(omni.time,omni.AE,superMag.time);

%%
% Selecting substorms closest to PFISR location
% deltaMLT = absDiffMLT(superMag.pfisrMlt,superMag.mlt);
% desiredMLTIndx = abs(deltaMLT)<Dmlt;
% desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat;
% closestSubstormIndx = desiredMLTIndx & desiredMLATIndx; 

%% Adding the PFISR experiments running during the substorm time

expBCArray = barker_coded_experiments();
% Barker Coded PFISR Experiment Filter
bcFilterIndx = zeros(1,length(amisr.expId));
for iexp = 1:1:length(expBCArray)
    bcFilterIndx = bcFilterIndx|strcmp(strtrim(expBCArray(iexp)),cellstr(deblank(amisr.expName)));
end

for iStorm = 1:1:length(superMag.stormID)
    tempIndx = find_amisr_exp(superMag.time(iStorm),amisr.startTime, amisr.endTime);
    numExp(iStorm) = length(tempIndx);
    amisrIndx(iStorm)=tempIndx(1);
end
superMag.expID = repmat(string("nan"),1,length(superMag.stormID));
superMag.expName = repmat(string("nan"),1,length(superMag.stormID));
superMag.status = repmat(string("nan"),1,length(superMag.stormID));
superMag.startTime = nan(1,length(superMag.stormID));
superMag.endTime = nan(1,length(superMag.stormID));
superMag.expBC = false(1,length(superMag.stormID));
superMag.numberOfSimultaneousExp = nan(1,length(superMag.stormID));

superMag.expID(~isnan(amisrIndx))=amisr.expId(amisrIndx(~isnan(amisrIndx)));
superMag.expName(~isnan(amisrIndx))=amisr.expName(amisrIndx(~isnan(amisrIndx)));
superMag.status(~isnan(amisrIndx))=amisr.status(amisrIndx(~isnan(amisrIndx)));
superMag.startTime(~isnan(amisrIndx))=amisr.startTime(amisrIndx(~isnan(amisrIndx)));
superMag.endTime(~isnan(amisrIndx))=amisr.endTime(amisrIndx(~isnan(amisrIndx)));
superMag.expBC(~isnan(amisrIndx))=bcFilterIndx(amisrIndx(~isnan(amisrIndx)));
superMag.numberOfSimultaneousExp(~isnan(amisrIndx))=numExp(amisrIndx(~isnan(amisrIndx)));


%% Create a table
T = table(superMag.datetime,superMag.stormID,...
    superMag.AE, superMag.mlat,...
    superMag.mlt, superMag.pfisrMlt,...
    superMag.expID',superMag.expName',...
    datetime(superMag.startTime,'ConvertFrom','datenum')',...
    datetime(superMag.endTime,'ConvertFrom','datenum')',...
    superMag.status',...
    superMag.expBC',...
    'VariableNames',{'Time','stormID','AE','MLAT','MLT','PFISR_MLT',...
    'PFISR_ExpID','PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
    'PFISR_ExpStatus','BarkerCode'});
%%
% disp('Table of Barker Code Experiments during a SuperMag substorm');
% T(T.BarkerCode,:)
%%
% disp('Table of PFISR Experiments during a SuperMag substorm');
% T(~strcmp(T.PFISR_ExpID,"nan"),:)

%% Finding if DASC is ON during a particular substorm with PFISR ON
% substormIndx = 1:1:length(superMag.time);
% cIndx = substormIndx(closestSubstormIndx);
% Tdasc = read_h5_data(dascFileStr);
% [timeStamp, wavelength] = restructure_DASC_table_to_time_array(Tdasc);
% 
% timeStamp = datenum(timeStamp);
% superMagCloseStorm.time = datetime(superMag.time(cIndx),'ConvertFrom','datenum');
% tempStart = datenum(dateshift(superMagCloseStorm.time,'start','day'));
% tempEnd = datenum(dateshift(superMagCloseStorm.time,'end','day'));
% superMagCloseStorm.DASC_timeMin = repmat(NaT,length(cIndx),1);
% superMagCloseStorm.DASC_timeMax = repmat(NaT,length(cIndx),1);
% superMagCloseStorm.DASC_wavelength = repmat({"nan"},length(cIndx),1);
% 
% Ft = griddedInterpolant(timeStamp,1:numel(timeStamp),'nearest','nearest');
% superMagCloseStorm.DASC_timeMin = datetime(timeStamp(Ft(tempStart)),'ConvertFrom','datenum');
% superMagCloseStorm.DASC_timeMax = datetime(timeStamp(Ft(tempEnd)),'ConvertFrom','datenum');
% tmIndex = superMagCloseStorm.DASC_timeMin==superMagCloseStorm.DASC_timeMax;
% superMagCloseStorm.DASC_timeMin(tmIndex) = NaT;
% superMagCloseStorm.DASC_timeMax(tmIndex) = NaT;
% 
% for cStorm = 1:1:length(cIndx)
%     if Ft(tempEnd(cStorm)) ~= Ft(tempStart(cStorm))
%     superMagCloseStorm.DASC_wavelength(cStorm,1) = {unique(wavelength(Ft(tempStart(cStorm)):Ft(tempEnd(cStorm))))};
%     end
% end
% 
%  %% Substorms that are close to PFISR
%  T1 = table(superMag.datetime(closestSubstormIndx),superMag.stormID(closestSubstormIndx),...
%     superMag.AE(closestSubstormIndx), superMag.mlat(closestSubstormIndx),...
%     superMag.mlt(closestSubstormIndx), superMag.pfisrMlt(closestSubstormIndx),...
%     superMag.expID(closestSubstormIndx)',superMag.expName(closestSubstormIndx)',...
%     datetime(superMag.startTime(closestSubstormIndx)','ConvertFrom','datenum'),...
%     datetime(superMag.endTime(closestSubstormIndx)','ConvertFrom','datenum'),...
%     superMag.status(closestSubstormIndx)',...
%     superMag.expBC(closestSubstormIndx)',...
%     superMagCloseStorm.DASC_timeMin,...
%     superMagCloseStorm.DASC_timeMax,...
%     superMagCloseStorm.DASC_wavelength,...
%     'VariableNames',...
%     {'Time','stormID','AE','MLAT','MLT','PFISR_MLT','PFISR_ExpID',...
%     'PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
%     'PFISR_ExpStatus','BarkerCode',...
%     'DASC_TimeMin','DASC_TimeMax','DASC_Wavelength'});

end

function [time, wavelength, url, wavelengthStr] = get_DASC_times_during_substorm(...
    dascTimeStamps, wavelengths, substormTime, growthDuration, expansionDuration)
    
    %substormTime - datetime
    if nargin<5
        expansionDuration = 1.0; %hr
    end
    
    if nargin<4
        growthDuration = 2.0; %hr
    end
    
    tempStart = substormTime - hours(growthDuration);
    tempEnd = substormTime + hours(expansionDuration);
    indxtStamp = dascTimeStamps<=tempEnd & dascTimeStamps>=tempStart;
    time = dascTimeStamps(indxtStamp);
    wavelength = wavelengths(indxtStamp);
    wavelengthStr = num2str(wavelength,'%04.f');
    url = create_DASC_url(time,wavelength);
end

function data=load_ascii_files(loadFile, format, headerlines)
if nargin<3
    headerlines = 1;
end
fileID = fopen(loadFile,'r');
data = textscan(fileID, format, 'headerlines', headerlines);
fclose(fileID);
end

function superMag = extract_superMag_data(loadFile, ftype)
if ftype == 1
    format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
    tempData = load_ascii_files(loadFile, format, 69);
    superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
    superMag.time = datenum(datestr(superMag.datetime));
    superMag.mlat = tempData{6};
    superMag.mlt = tempData{7};
end
end

function amisr = extract_amisr_data(loadFile)
    table = read_h5_data(loadFile);
    amisr.expId = string(table.Data{2});
    amisr.expName = string(table.Data{3});
    amisr.status = string(table.Data{4});
    amisr.startTime = unixtime2matlab(table.Data{5});
    amisr.endTime = unixtime2matlab(table.Data{1});
end

function [Lm,MLT] = get_pfisr_magnetic_coordinates(time,maginput,GDZ,magFieldNo)
    [Lm,~,~,~,~,MLT] = onera_desp_lib_make_lstar(magFieldNo,[0,0,0,0,0],0,time,GDZ(3),GDZ(1),GDZ(2),maginput);
%     tup=py.aacgmv2.wrapper.get_aacgm_coord(GDZ(1), GDZ(2), GDZ(3), time, 'TRACE');
%     MLAT = double(py.array.array('d',py.numpy.nditer(tup{1})));
%     MLON = double(py.array.array('d',py.numpy.nditer(tup{2})));
%     MLT_AACGM = double(py.array.array('d',py.numpy.nditer(tup{3})));
    Lm = abs(Lm);
end

function [amisrIndx] = find_amisr_exp(time, startTimeArr, endTimeArr)
% Finds the amisr experiment indx
    amisrIndx = find(time>startTimeArr & time<endTimeArr);
    if isempty(amisrIndx)
        amisrIndx=nan;
    end
    
end

function expArr = barker_coded_experiments()
expArr =["GenPINOT_PulsatingAurora_TN30          ";
    "Inspire_v01                            ";
    "Kelley01                               ";
    "MSWinds23                              ";
    "MSWinds23_dt013                        ";
    "MSWinds23_3dt                            ";
    "MSWinds23m                               ";
    "MSWinds23hr                              ";
    "MSWinds21                                ";
    "MSWinds26.v03                          ";
    "Semeter01                              ";
    "Sporadic01                             ";
    "Sporadic02                             ";
    "Sporadic03                             ";
    "Sporadic04                             ";
    "Sporadic14                             ";
    "Sporadic15                             ";
    "Sporadic15_3dt                         ";
    "ThemisD1.v01                             "
    "MSWinds27.v01                            "];
end
 
function [timeStamp, wavelength] = restructure_DASC_table_to_time_array(T)
    wavelength = string();
    timeStamp=cell2mat(T.Data(strcmp(string(T.Name),"time")));
    waveCode=T.Data(strcmp(string(T.Name),"wavelengthCode"));
    waveLength=T.Data(strcmp(string(T.Name),"wavelength"));
    for i = 0:1:length(waveCode)-1
        if i==0
            wavelength(end:end+numel(waveLength{i+1})-1,1) = string(waveCode{i+1}(waveLength{i+1},1));
        else
            wavelength(end+1:end+numel(waveLength{i+1}),1) = string(waveCode{i+1}(waveLength{i+1},1));
        end
    end
    [timeStamp,I] = sort(timeStamp);
    [timeStamp, w] = unique(timeStamp, 'stable');
    wavelength = wavelength(I);
    wavelength = wavelength(w);
    timeStamp = datetime(unix_to_matlab_time(timeStamp),'ConvertFrom','datenum');
end

function dMLT = absDiffMLT(a,b)
    % Source: https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles 
    normMLT = mod(a-b,24);
    dMLT = min(24-normMLT, normMLT);
end
