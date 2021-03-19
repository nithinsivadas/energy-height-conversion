clear all; 

storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
load([storeDir,'table_of_substorms_as_input.mat']);
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
workDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir';

%%
fileNameList = struct2cell(dir([workDir,strcat(filesep,'*_pfisrData.h5')]));
filePathStr = strcat(strcat(workDir,filesep),string(fileNameList(1,:)'));

for i=1:1:length(filePathStr)
    storedStormID(:,i) = str2double(fileNameList{1,i}(1:regexp(fileNameList{1,i},'_')-1));
end

%%
time = unixtime2matlab(h5read(omniFile,'/Time'));
SML = h5read(omniFile,'/Indices/SML');
SMU = h5read(omniFile,'/Indices/SMU');
Pdyn = h5read(omniFile,'/FlowPressure');
Fp = griddedInterpolant(time,Pdyn);
Fsml = griddedInterpolant(time,SML);
Fsmu = griddedInterpolant(time,SMU);
%%
BzGSM = h5read(omniFile,'/BField/BzGSM');
FBz = griddedInterpolant(time,BzGSM);
%%
ByGSM = h5read(omniFile,'/BField/ByGSM');
FBy = griddedInterpolant(time,ByGSM);

%% Magnetopause parameters from Shen et al., 1993

r_0 = (11.4 + 0.013.*BzGSM).*(Pdyn).^(-1./6.6); % for Bz>=0
r_0(BzGSM<0) = (11.4 + 0.14.*BzGSM(BzGSM<0)).*(Pdyn(BzGSM<0)).^(-1./6.6); % for Bz<0
% r_0 is the standoff distance in R_E
% alpha is the level of tail flaring
alpha = (0.58-0.010*BzGSM).*(1+0.01*Pdyn);
Fr_0 = griddedInterpolant(time,r_0);
Falpha = griddedInterpolant(time,alpha);

%%
BxGSE = h5read(omniFile,'/BField/BxGSE');
ByGSE = h5read(omniFile,'/BField/ByGSE');
BzGSE = h5read(omniFile,'/BField/BzGSE');
B = (BxGSE.^2+ByGSE.^2+BzGSE.^2).^0.5;
FB = griddedInterpolant(time,B);
%%
AL = h5read(omniFile,'/Indices/AL');
AL(AL==99999)=nan;
FAL = griddedInterpolant(time, AL);
%%
density = h5read(omniFile,'/ProtonDensity');
Fdensity = griddedInterpolant(time,density);
%%
velocity = h5read(omniFile,'/Velocity/V');
Fv = griddedInterpolant(time,velocity);

%%
MLT_qualification = T.MLT>=16 | T.MLT<=4;
T = T(MLT_qualification,:);
isolatedSubstormDuration = repmat(duration,length(T.Time),1);
isolatedSubstormDuration(2:end) = T.Time(2:end) - T.Time(1:end-1);

%%
BzGSM(BzGSM==0) = 0.0001;
B_T = (ByGSM.^2 + BzGSM.^2).^0.5;

theta_c = wrapTo2Pi(atan2(ByGSM,BzGSM));
Ftheta = griddedInterpolant(time,theta_c);
%%
l_0 = 7*(6371*10^3);
E_kl = velocity.*B_T.*(sin(theta_c/2)).^2; %Km nT/s

E = 1e-9.*1e7.*(velocity.*10.^3).*((B.*10^-9).^2).*(sin(theta_c/2)).^4*l_0^2; %GW

Fekl = griddedInterpolant(time,E_kl);
FE = griddedInterpolant(time,E);
%%
T.previousSubstormDuration = isolatedSubstormDuration;
T.nextSubstormDuration = repmat(duration,length(T.Time),1);
T.nextSubstormDuration(1:end-1) = T.Time(2:end) - T.Time(1:end-1);
%%
for i=1:1:length(T.stormID)
      
    if sum(storedStormID==T.stormID(i))>0
        T.storageLocation(i) = filePathStr(storedStormID==T.stormID(i));
    end
    
    tempTime = T.Time(i)-duration(1,0,0):duration(0,1,0):T.Time(i)+duration(0,30,0);
    [T.SML(i),indx2] = min(Fsml(datenum(tempTime)));
    T.peakSMLTime(i)= tempTime(indx2);
    
end

%%
timeMinStr = "01 Dec 2006";
timeMaxStr = "31 Dec 2019";


% Print out substorms that might have the radiation belt signatures

%%
tic
% T5=T4(T4.previousSubstormDuration>duration(3,0,0) & T4.MLAT<65.8 ,:); %& ~ismissing(T4.storageLocation)
% T5=T1(T1.previousSubstormDuration>duration(3,0,0),:);
T5 = T;
% T5=T4;
GP = duration(3,0,0);
EP = duration(3,0,0);
smlTime = -GP:duration(0,1,0):EP;
T5.timeArray = repmat(duration,height(T5),length(smlTime));
T5.smlArray = zeros(height(T5),length(smlTime));
T5.smuArray = zeros(height(T5),length(smlTime));
T5.ALArray = zeros(height(T5),length(smlTime));
T5.BzArray = zeros(height(T5),length(smlTime));
T5.pArray = zeros(height(T5),length(smlTime));
T5.densityArray = zeros(height(T5),length(smlTime));
T5.vArray = zeros(height(T5),length(smlTime));
T5.EklArray = zeros(height(T5),length(smlTime));
T5.EArray = zeros(height(T5),length(smlTime));
T5.thetaArray = zeros(height(T5),length(smlTime));
T5.r_0Array = zeros(height(T5),length(smlTime));
T5.alphaArray = zeros(height(T5),length(smlTime));

for i=1:1:height(T5)
    TMatrix(i,:) = datenum(T5.Time(i)-GP:duration(0,1,0):T5.Time(i)+EP);
end
T5.timeArrayRelOnset = repmat(smlTime,height(T5),1);
T5.timeArray = TMatrix;

T5.smlArray = Fsml(TMatrix);
T5.smuArray = Fsmu(TMatrix);
T5.pArray   = Fp(TMatrix); 
T5.BzArray  = FBz(TMatrix); 
T5.ALArray  = FAL(TMatrix); 
T5.densityArray = Fdensity(TMatrix); 
T5.vArray   = Fv(TMatrix); 
T5.eklArray = Fekl(TMatrix); 
T5.EArray   = FE(TMatrix); 
T5.thetaArray   = Ftheta(TMatrix); 
T5.r_0Array   = Fr_0(TMatrix); 
T5.alphaArray = Falpha(TMatrix);

toc

%%
condition = T5.previousSubstormDuration>duration(6,0,0)... %    
    & T5.nextSubstormDuration>duration(3,0,0)...
    & ~ismissing(T5.storageLocation)...
    & absDiffMLT(T5.MLT, T5.PFISR_MLT)<2 ...
    & T5.MLAT<66;

T6 = T5(condition,:);

altD = 60:1:97;
alt = 60:1:140;
for i=1:1:size(T6,1)
    
    pfisrNe = h5read(T6{i,'storageLocation'},'/inputData/Ne');
    pfisrHall = h5read(T6{i,'storageLocation'},'/conductivity/hall');
    pfisrEflux = h5read(T6{i,'storageLocation'},'/energy/energyFlux');
    pfisrEbins = h5read(T6{i,'storageLocation'},'/energy/energyBin'); 
    pfisrtime = h5read(T6{i,'storageLocation'},'/time');
    pfisralt =  h5read(T6{i,'storageLocation'},'/alt');
    pfisrRelTime = interp1(T6.timeArray(i,:),T6.timeArrayRelOnset(i,:),pfisrtime,'linear','extrap');
    T6.pfisrNe{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrNe,alt,'linear')',T6.timeArrayRelOnset(i,:));
    T6.pfisrHall{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrHall,alt,'linear')',T6.timeArrayRelOnset(i,:));
    T6.pfisrEbins{i} = pfisrEbins;
    T6.pfisrEflux{i} = interp1(pfisrRelTime,pfisrEflux',T6.timeArrayRelOnset(i,:));

    
end
%%

h=figure; 
p=create_panels(h,'totalPanelNo',10,'marginbottom',10,'panelHeight',40);

% condition = T5.previousSubstormDuration>duration(3,0,0)...
%     & T5.nextSubstormDuration>duration(3,0,0)...
%     & T5.SML<-800;
%     & ~ismissing(T5.storageLocation)...
%     & absDiffMLT(T5.MLT, T5.PFISR_MLT)<2 ...
%     & T5.MLAT<66;
%
% condition1 = T5.previousSubstormDuration>duration(3,0,0)...
%     & T5.nextSubstormDuration>duration(3,0,0)...
%     & T5.SML<-800;
% 
% condition2 = T5.previousSubstormDuration>duration(3,0,0)...
%     & T5.nextSubstormDuration<duration(3,0,0)...
%     & T5.SML<-800;
%
p(1,1).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'BzArray'});
ylabel('IMF Bz [nT]')
ylim([-3 1]);
set(gca,'XTick',{});
hold on;
plot3([min(smlTime) max(smlTime)],[nanmedian(BzGSM) nanmedian(BzGSM)],[1 1]);
hold on;
plot3([duration duration],[-3 1],[1 1]);
% caxis([0,0.4]);
title(['No. of Substorms: ',num2str(length(T5.stormID(condition)))]);

p(1,2).select();
time_plot(T5.timeArrayRelOnset(1,:),sin(T5{condition,'thetaArray'}./2).^2);
ylabel('sin^2(\Theta/2)');
ylim([0.4 0.8]);
set(gca,'XTick',{});
hold on;
plot3([min(smlTime) max(smlTime)],[nanmedian((sin(theta_c./2)).^2) nanmedian((sin(theta_c./2)).^2)],[1 1]);
hold on;
plot3([duration duration],[0 1],[1 1]);
% caxis([0,4]);

p(1,3).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'pArray'});
ylabel('P_{dyn} [nPa]');
ylim([1.6,2.8]);
set(gca,'XTick',{});
hold on;
plot3([min(smlTime) max(smlTime)],[nanmean(Pdyn) nanmean(Pdyn)],[1 1]);
hold on;
plot3([duration duration],[1.6 2.8],[1 1]);
% caxis([0,1.5]);


p(1,4).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'smlArray'});
ylabel('SML [nT]')
ylim([-400,0]);
set(gca,'XTick',{});
hold on;
plot3([min(smlTime) max(smlTime)],[nanmedian(SML) nanmedian(SML)],[1 1]);
hold on;
plot3([duration duration],[-400 0],[1 1]);
% caxis([0 0.015]);

p(1,5).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'r_0Array'});
ylabel('r_0 [R_E]')
ylim([9.5,10.5]);
set(gca,'XTick',{});
hold on;
plot3([min(smlTime) max(smlTime)],[nanmedian(r_0) nanmedian(r_0)],[1 1]);
hold on;
plot3([duration duration],[9 11],[1 1]);
% caxis([0 1.5]);

p(1,6).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'alphaArray'});
ylabel('\alpha [a.u.]')
% ylim([9.5,10.5]);
set(gca,'XTick',{});
% hold on;
% plot3([min(smlTime) max(smlTime)],[nanmedian(r_0) nanmedian(r_0)],[1 1]);
% hold on;
% plot3([duration duration],[9 11],[1 1]);
% caxis([0 1.5]);

p(1,7).select();
time_plot(T5.timeArrayRelOnset(1,:),T5{condition,'EArray'});
ylabel('\epsilon(t) [GW]');

hold on;
plot3([min(smlTime) max(smlTime)],[nanmedian(E) nanmedian(E)],[1 1]);
hold on;
plot3([duration duration],[0 600],[1 1]);
% set(gca,'YScale','log');
ylim([10 400]);
set(gca,'XTick',{});



p(1,8).select();
time_plot_2D(T6.timeArrayRelOnset(1,:),T6{:,'pfisrNe'},alt);
set(gca,'ColorScale','log','CLim',[10^10 10^12]);
set(gca,'XTick',{});
ylabel('N_e');

p(1,9).select();
time_plot_2D(T6.timeArrayRelOnset(1,:),T6{:,'pfisrEflux'},T6{:,'pfisrEbins'}{1});
set(gca,'ColorScale','log','CLim',[10^8 10^12],'YScale','log');
xlim([-3,3]);
set(gca,'XTick',{});
ylabel('\phi(E)');

p(1,10).select();
time_plot_2D(T6.timeArrayRelOnset(1,:),T6{:,'pfisrHall'},alt);
set(gca,'ColorScale','log','CLim',[10^-6 10^-2]);
ylabel('\Sigma_H');
xlabel('Substorm Time (t)');
% caxis([0,0.01]);
% figure; 
% time_plot(smlTime,((vArray/1000).^2.*(densityArray*10^6)).^(-1/6));
% ylabel('magnetopause distance [a.u.]');
% ylim([0.09,0.11]);

% figure; 
% time_plot(smlTime,AlArray);
% ylabel('AL [nT]');
% 
%%
% figure; 
% time_plot(smlTime,densityArray);
% ylabel('Density [cm^{-3}]');
% ylim([0,10]);
%%
% figure; 
% time_plot(smlTime,pArray);
% ylabel('Dynamic Pressure');
% ylim([1.6,2.8]);
%%

% figure; 
% time_plot(smlTime,vArray);
% ylabel('Velocity [km/s]');
% ylim([350 650]);
% 

%%
% figure; 
% time_plot_2D(T6.timeArrayRelOnset(1,:),T6{:,'pfisrNe'},alt);
% set(gca,'ColorScale','log','CLim',[10^9 10^12]);
% datetick('x','HH:MM');
% xticks(T6.timeArrayRelOnset(1,:));
% ylabel('IMF Bz [nT]');
% ylim([-2 0.5]);

%%

%%

%%
function time_plot(t,y)


% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;
my = nanmean(y);
medy = nanmedian(y);
sdev = nanstd(y);
sy = sdev./((size(y,1).^0.5));

edges = linspace(min(my)-nanmean(sdev),max(my)+nanmean(sdev),300);
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

function dMLT = absDiffMLT(a,b)
    % Source: https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles 
    normMLT = mod(a-b,24);
    dMLT = min(24-normMLT, normMLT);
end

