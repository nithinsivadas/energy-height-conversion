%% Finding why E(V|I) is different from E(I|V)
clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%
omniFile = 'G:\My Drive\Research\Projects\Data\omni_test.h5';
omni = extract_omni_data(omniFile);
%%
C = define_universal_constants;
Time = datenum('01 Jan 1981'):1/1440:datenum('01 Jan 2019');
E = omni.Fekl(Time).*10^-3;
J = (2./C.mu0).*omni.Fsml(Time).*10^-6;
% figure; 
% 
% scatter(E,-J,10,'filled','MarkerFaceAlpha',0.2);
% xlim([0,40]);
% ylim([0,6000]);
% xlabel(' E_M [mV/m]');
% ylabel('J_{iono}(SML) [mA/m^2]]');

E_signal = movmean(E,20);
E_noise = movvar(E,20);
E_snr = E_signal./E_noise; 

Ebins = 0:1:40;
for i = 1:1:length(Ebins)-1
%     ESNR(i) = nanmean(E_snr(E>=Ebins(i) & E<Ebins(i+1)));
    ESIGNAL(i) = nanmean(E(E_signal>=Ebins(i) & E_signal<Ebins(i+1)));
    ENOISE(i) = nansum(E_noise(E_signal>=Ebins(i) & E_signal<Ebins(i+1)));
end

figure; 
plot(Ebins(1:end-1),ESIGNAL);
ylabel('E_{SIGNAL}');

figure; 
plot(Ebins(1:end-1),ENOISE);
ylabel('E_{NOISE}');

ESNR = ESIGNAL./ENOISE;
figure; 
plot(Ebins(1:end-1),ESNR);
ylabel('E_{SNR}');


J_signal = movmean(-J,180);
J_noise = movvar(-J,180);
J_snr = J_signal./J_noise; 

Jbins = linspace(0,6000,100);
for i = 1:1:length(Jbins)-1
%     JSNR(i) = nanmean(J_snr(-J>=Jbins(i) & -J<Jbins(i+1)));
    JSIGNAL(i) = nanmean(-J(J_signal>=Jbins(i) & J_signal<Jbins(i+1)));
    JNOISE(i) = nansum(J_noise(J_signal>=Jbins(i) & J_signal<Jbins(i+1)));
%     JNOISE(i) = nanvar(-J(J_signal>=Jbins(i) & J_signal<Jbins(i+1)));
end

JSNR = JSIGNAL./JNOISE;
figure; 
plot(Jbins(1:end-1),JSIGNAL./JNOISE);
ylabel('J_{SNR}');

figure;
plot(Jbins(1:end-1),JSIGNAL);
ylabel('J_{Signal}');

figure;
plot(Jbins(1:end-1),JNOISE);
ylabel('J_{Noise}');


%% Finding mean and variance of E and J

clearvars selection ESIGNAL ENOISE 

C = define_universal_constants;
Time = datenum('01 Jan 1981'):1/1440:datenum('01 Jan 2019');
E = omni.Fekl(Time).*10^-3;
J = (2./C.mu0).*omni.Fsml(Time).*10^-6;

tic
[E_signal, E_sample_matrix, index] = runmean(E,20);
toc

Ebins = linspace(0,40,100);
% Ebins = 0:1:40;


for i = 1:1:length(Ebins)-1
%     ESNR(i) = nanmean(E_snr(E>=Ebins(i) & E<Ebins(i+1)));
    selection(i,:) = E_signal>=Ebins(i) & E_signal<Ebins(i+1);
    ESIGNAL(i) = nanmean(E(selection(i,:)));
%     ENOISE(i) = nansum(E_noise(E_signal>=Ebins(i) & E_signal<Ebins(i+1)));
end

for i = 1:1:length(Ebins)-1
    
    ENOISE(i) = nanvar(reshape(E_sample_matrix(selection(i,:),:),[],1));
   
end

figure; 

plot(Ebins(1:end-1), ESIGNAL./ENOISE);
xlabel('E [mV/m]');
ylabel('mean(E)/var(E)');


%%
clearvars selection JSIGNAL JNOISE 

C = define_universal_constants;
Time = datenum('01 Jan 1981'):1/1440:datenum('01 Jan 2019');
E = omni.Fekl(Time).*10^-3;
J = (2./C.mu0).*omni.Fsml(Time).*10^-6;

tic
[J_signal, J_sample_matrix, index] = runmean(-J,20);
toc

Jbins = linspace(0,6000,100);


for i = 1:1:length(Jbins)-1
%     ESNR(i) = nanmean(E_snr(E>=Ebins(i) & E<Ebins(i+1)));
    selection(i,:) = J_signal>=Jbins(i) & J_signal<Jbins(i+1);
    JSIGNAL(i) = nanmean(-J(selection(i,:)));
%     ENOISE(i) = nansum(E_noise(E_signal>=Ebins(i) & E_signal<Ebins(i+1)));
end

for i = 1:1:length(Jbins)-1
    
    JNOISE(i) = nanvar(reshape(J_sample_matrix(selection(i,:),:),[],1));
   
end

%%
figure; 

plot(Jbins(1:end-1), JSIGNAL./JNOISE);
xlabel('J [mV/m]');
ylabel('mean(J)/var(J)');

% Case 1: noisy V and I data (gaussian)
%%
clearvars Ns Ns2 EIgV EVgI
SNR = 10; % In linear unitsl and is equal to (rms(Vs)./rms(Ns)).^2
SNRArray = logspace(0.5,-1.5,1000).*SNR;
SNRArray1 = interp1(Ebins(1:end-1),ESNR,linspace(Ebins(1),Ebins(end-1),1000));
SNRArray2 = interp1(Jbins(1:end-1),JSNR,linspace(Jbins(1),Jbins(end-1),1000));

signals = linspace(0,20,1000); % DC signal voltage for 500 signals
nRepeat = 1000;
Vs = repmat(signals,nRepeat,1); %500 instances of each signal
for i = 1:length(signals)
%     Ns(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(0.01.*SNRArray(i));
    Ns(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(SNRArray1(i));
end

for i = 1:length(signals)
%     Ns2(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(0.03.*SNRArray(i));
    Ns2(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(SNRArray2(i));
end

yLimRange = [0,6000];
xLimRange = [0,40];

% figure; plot(rms(Vs).^2./rms(Ns).^2);
V = Vs + Ns;
% V = Vs;

% R = linspace(20,1,500); %Constant
R = 1/300;

Isignal = (signals./R) ;
Is = (Vs./R);
INs = (Ns./R);
% I = (V./R);
I = Is + Ns2./R;
% I = Is + INs;

vRange = [0,max(V(:))];
iRange = [0,max(I(:))];

vBinsTemp = linspace(0,max(V(:)),2000);
vBins = vBinsTemp(1:end-1) + diff(vBinsTemp);

for i = 1:1:length(vBinsTemp)-1
    EIgV(i) = mean(I(V>=vBinsTemp(i) & V<vBinsTemp(i+1) & V >=0 )); 
end


iBinsTemp = linspace(0,max(I(:)),2000);
iBins = iBinsTemp(1:end-1) + diff(iBinsTemp);

for i = 1:1:length(iBinsTemp)-1
    EVgI(i) = mean(V(I>=iBinsTemp(i) & I<iBinsTemp(i+1) & I >=0 )); 
end

hFig = figure; 
resize_figure(hFig,80,250);

    p = panel(hFig);
    p.pack({{60}},3);
    p.marginright=20;

    p(1,1).select();
    scatter(V(:),I(:),10,'filled','MarkerFaceAlpha',0.01);
    r = corrcoef(V(:),I(:));
    text(0.8,0.1,['r_0_0 = ',num2str(r(2,1),3)],'Units','normalized');
    xlabel('E');
    ylabel('J');
    hold on;
    plot(signals,Isignal);
    xlim(xLimRange);
    ylim(yLimRange);
    legend('Signal + Noise','Signal');
    
    p(1,2).select();
    yyaxis left
    plot(vBins,EIgV); 
    ylabel('\langle J|E \rangle');
    xlim(xLimRange);
    ylim(yLimRange);
%     xlim(vRange);
%     ylim(iRange);
    hold on; 

    yyaxis right
    plot(signals,Isignal);
    ylabel('J');
%     ylim(iRange);
    ylim(yLimRange);
    xlabel('E');

    legend('Signal + Noise', 'Signal');

    p(1,3).select(); 
    plot(EVgI,iBins); 
    hold on; 
    plot(signals,Isignal);
    xlim(xLimRange);
    ylim(yLimRange);
%     xlim(vRange);
%     ylim(iRange);
    xlabel('\langle E|J \rangle');
    ylabel('J');
    legend('Signal + Noise', 'Signal','Location','northeast');
    
    
%% Calculating ACE and WIND data and their errors 

ace = extract_data('G:\My Drive\Research\Projects\Data\omni_components\ace_min_b2017.txt');
dscov = extract_data('G:\My Drive\Research\Projects\Data\omni_components\dscov_min_b2017.txt');
wind = extract_data('G:\My Drive\Research\Projects\Data\omni_components\wind_min_b2017.txt');

%%
figure;
[h,xedge,yedge]=histcounts2((ace.E_kl-wind.E_kl)*10^-3,wind.E_kl*10^-3,-10:0.1:10,0:0.3:20,'Normalization','pdf');
[h1,xedge1]=histcounts(wind.E_kl*10^-3,0:0.3:20,'Normalization','pdf');
[X,Y] = meshgrid(xedge(1:end-1),yedge(1:end-1));
% figure; 
p=pcolor(X,Y,(h./h1)');
% p=pcolor(X,Y,h');
% caxis([0,0.3]);
set(p, 'EdgeColor', 'none');
colorbar;
% % figure;
% scatter((ace.E_kl-wind.E_kl)*10^-3,wind.E_kl*10^-3,'filled');

%% Estimating kernal density

gride = -5:0.1:5;
gridx = 0:0.5:20;
[E, X] = meshgrid(gride,gridx);

[fx,xi] = ksdensity((wind.E_kl*10^-3),gridx);
Fx = griddedInterpolant(xi,fx);
[fex,xii] = ksdensity([(ace.E_kl-wind.E_kl)*10^-3, ace.E_kl*10^-3],[E(:),X(:)]);
Fex = scatteredInterpolant(xii(:,1),xii(:,2),fex);

% gride = -40:1:40;
% gridx = 300:1:800;
% [E, X] = meshgrid(gride,gridx);

% [fx,xi] = ksdensity((wind.velocity),gridx);
% Fx = griddedInterpolant(xi,fx);
% [fex,xii] = ksdensity([(ace.velocity-wind.velocity), wind.velocity],[E(:),X(:)]);
% Fex = scatteredInterpolant(xii(:,1),xii(:,2),fex);
Dx = diff(gridx);
Dx(end+1) = Dx(end);
%%

gride1 = -5:0.01:5;
gridx1 = 0:0.1:20;
[E1, X1] = meshgrid(gride1,gridx1);
figure; 
p=pcolor(E1,X1,(Fex(E1,X1)./(Fx(gridx1)')./sum((Fex(E1,X1)./(Fx(gridx1)')),2)));
set(p, 'EdgeColor', 'none');
colorbar;
caxis([0,0.05]);
%%

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

  function [y,D,index] = runmean(x,w)
    if rem(w,2)
        error('Only works for even window values');
    end
    
    m = round(0.5*(w-1)); 
    k=1;
    y = zeros(length(1+m:length(x)-m),1);
    D = zeros(length(1+m:length(x)-m),w);
    index = 1+m:length(x)-m;
    for i = index
        j=i-m:i+m-1;
        y(k) = nansum(x(j))./w;
        D(k,:) = x(j);
        k=k+1;
    end
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
    
    symH = h5read(omniFile,'/Indices/SYM_H');
    omni.FsymH = griddedInterpolant(time, symH);
    % Solar wind
    
%         % Dynamic Pressure
%     Pdyn = h5read(omniFile,'/FlowPressure');
%     omni.Fp = griddedInterpolant(time,Pdyn);
    
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
    
    EField = h5read(omniFile,'/EField');
    omni.FEfield = griddedInterpolant(time,EField);
    
    E = 1e-9.*1e7.*(velocity.*10.^3).*((B.*10^-9).^2).*(sin(theta_c/2)).^4*l_0^2; %GW 
    omni.FE = griddedInterpolant(time,E);
    
    
%     % Magnetopause parameters from Shen et al., 1993
% 
%         % r_0 is the standoff distance in R_E
%     r_0 = (11.4 + 0.013.*BzGSM).*(Pdyn).^(-1./6.6); % for Bz>=0
%     r_0(BzGSM<0) = (11.4 + 0.14.*BzGSM(BzGSM<0)).*(Pdyn(BzGSM<0)).^(-1./6.6); % for Bz<0
%     omni.Fr_0 = griddedInterpolant(time,r_0);
%     
%         % alpha is the level of tail flaring
%     alpha = (0.58-0.010*BzGSM).*(1+0.01*Pdyn);
%     omni.Falpha = griddedInterpolant(time,alpha);
    
    
end


