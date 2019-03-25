%% Contribution of optical emissions from precipitation (GLOW)
% Paper 2
clear all
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';

%% Load data
pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');

pfisrImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);

pfisrImageNumF = pfisrImage./permute(repmat(pfisrData.zEnergyBin,1,1,size(pfisrImage,1)),[3,1,2]); % m-2 sr s eV
%% 
timeMinStr = '26 Mar 2008 11:04';
timeMaxStr = '26 Mar 2008 11:25';

[data,time,alt] = get_glow_parameter(timeMinStr, timeMaxStr, pfisrData, pfisrImageNumF,12,1e3,1e6);

%%
figure;
plot_2D_time_series(time,alt,data.A4278',5/60,0);
label_time_axis(time,true,5/60);
% ylim([60,150]);
colorbar;
%%
figure;
iB = 12;
flux = 2*pi*squeeze(pfisrImageNumF(:,iB,:))'/(1e4);
flux(flux<=0) = 0.000001;
plot_2D_time_series(pfisrData.time,pfisrData.zEnergyBin(iB,:),log10(flux),5/60,3,'26 Mar 2008 11:00','26 Mar 2008 11:30');
label_time_axis(time,true,5/60,'26 Mar 2008 11:00','26 Mar 2008 11:30');
caxis([1,4]);
colorbar;

%%
figure; 
time = datenum(2008,03,26,11,30,0);
glat = 65.1;
glon = -147.5;
[f107a, f107, aph] = f107_aph(time);
Ap = aph(:,1);
f107p = (f107a + f107)./2;
Echar = 100e3;
% Energy Grid
Emin = 1;
Emax = 1e6;
Nbins = 250;
[Ebins, Phitop] = monoenergetic(Emin, Emax, Nbins, Echar); 
% Phitop = energetic electron differential number flux into top of atmosphere: cm-2 s-1 eV-1
Phitop = Phitop*10;
iono0 = glowenergy(time, glat, glon, f107a, f107, f107p, Ap, Ebins, Phitop*0);
iono = glowenergy(time, glat, glon, f107a, f107, f107p, Ap, Ebins, Phitop);
%%
semilogx(iono.A5577-iono0.A5577,iono.altkm); %%cm-3 s; Rayleigh is 10^6 cm-2 s (integrated along altitude)
ylim([60 200]);

function [data,time,alt] = get_glow_parameter(timeMinStr, timeMaxStr, pfisrData, pfisrImageNumF, iB, Emin, Emax)

timeMinIndx = find_time(pfisrData.time, timeMinStr);
timeMaxIndx = find_time(pfisrData.time, timeMaxStr);
eMinIndx = find_altitude(pfisrData.zEnergyBin(iB,:),Emin);
eMaxIndx = find_altitude(pfisrData.zEnergyBin(iB,:),Emax);
nEnergy = length(pfisrData.zEnergyBin(1,:));
Ebins = logspace(log10(1),log10(pfisrData.zEnergyBin(1,1)),30);
Ebins(end+1:end+nEnergy-1) = pfisrData.zEnergyBin(1,2:end);
nBeams = length(pfisrData.zEnergyBin(:,1));
glat = mean(pfisrData.latitude(:));
glon = mean(pfisrData.longitude(:));
[f107a, f107, aph] = f107_aph(mean(pfisrData.time(timeMinIndx:timeMaxIndx)));
Ap = aph(:,1);
f107p = (f107a + f107)./2;

Phitop = zeros(size(Ebins));

k=1;
multiWaitbar('Generating optical emissions',0);
dt = 1./length(timeMinIndx:timeMaxIndx);
for iTime = timeMinIndx:timeMaxIndx
        multiWaitbar('Generating optical emissions','Increment',dt);
        flux = squeeze(pfisrImageNumF(iTime,iB,:))*(1e-4)*2*pi; %cm-2 s eV (Half solid angle)
        flux1 = zeros(size(flux));
        flux1(eMinIndx:eMaxIndx) = flux(eMinIndx:eMaxIndx);
        Phitop(end-nEnergy+1:end)=flux1;
        iono = glowenergy(pfisrData.time(iTime), glat, glon, f107a, f107, f107p, Ap, Ebins, Phitop);
        data.A5577(k,:) = iono.A5577;
        data.A4278(k,:) = iono.A4278;
        data.ionrate(k,:) = iono.ionrate;
        time(k) = pfisrData.time(iTime);
        k = k+1;
end
alt = iono.altkm;
end