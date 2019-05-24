%% Estimate the optical intensity contribution of PFISR estimated energy spectra
clear all;
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20101018.001_bc_15sec-full_v110.h5';
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
h5FileStrNe = 'G:\My Drive\Research\Projects\Paper 2\Data\Event 1\20080326.001_bc_15sec-fitcal.h5';
%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');
%%
pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');
%%
timeMinStr = '26 Mar 2008 10:30';
timeMaxStr = '26 Mar 2008 11:40';
%%
dascImage = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
pfisrEnergyFluxImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
pfisrNeImage = permute(h5read(h5FileStr,'/inputData/Ne'),[3 2 1]);
%%
pfisrNeAlt = h5read(h5FileStr,'/energyFluxFromMaxEnt/alt')';
%%
beamCodes = (h5read(h5FileStrNe,'/BeamCodes'));
%%
% magBeamNo = find_field_aligned_beam_no(beamCodes,-154.3,77.5);
magBeamNo = 12; % In order to match with Fig 4c; which picks latitude: 65.415N. Beam 12, is latitudinally close for 30 keV. 
%% Glow
time = pfisrData.time;
energyBin = pfisrData.zEnergyBin(magBeamNo,:);
numFlux = energy_to_num(squeeze(pfisrEnergyFluxImage(:,magBeamNo,:))',time,energyBin)'.*(10^-4); %cm-2 s-1 eV-1
numFlux(numFlux<0) = 0;
electronDensity = squeeze(pfisrNeImage(:,magBeamNo,:))'*10^-6;
%%
glat = pfisrData.latitude(magBeamNo,1);
glon = pfisrData.longitude(magBeamNo,1);
%%
nBins = 350; 
% energyBinGlow = logspace(1,6,nBins);
%%
energyBinGlow = egrid(10^6, 350);
%%
filter = ones(1,nBins);
filterLE = ones(1,nBins); %LE < 30 keV
filterLE(energyBinGlow>30000) = 0;
filterHE = ones(1,nBins); %HE >30 keV
filterHE(energyBinGlow<30000) = 0;
filterSE = ones(1,nBins); % SE - >100 keV
filterSE(energyBinGlow<100000) = 0;
% energyBinGlowLE = logspace(1,log10(30*10^3),250);
% energyBinGlowHE = logspace(log10(30*10^3),6,250);
%%
[f107a,f107,aph] = f107_aph(floor(time(1)));
Ap = aph(1);
tic
k=1;
multiWaitbar('Progress',0);
timeIndx = find_time(time,timeMinStr):1:find_time(time,timeMaxStr);
nTime = length(timeIndx);
dt = 1./nTime;
for iTime = timeIndx
      multiWaitbar('Progress','Increment',dt);
      iono = calculate_rayleigh(time(iTime), glat, glon, Ap, energyBin, numFlux(iTime,:), energyBinGlow, filter);
      glow.Ne(k,:) = iono.NeCalc;
      glow.NeIRI(k,:) = iono.Ne;
      glow.A5577(k,:) = iono.A5577;
      glow.A4278(k,:) = iono.A4278;
      glow.R5577(k,:) = iono.R5577;
      glow.R4278(k,:) = iono.R4278;    
      glow.A6300(k,:) = iono.A6300;
      glow.R6300(k,:) = iono.R6300;

      ionoLE = calculate_rayleigh(time(iTime), glat, glon, Ap, energyBin, numFlux(iTime,:), energyBinGlow, filterLE);
      glowLE.Ne(k,:) = ionoLE.NeCalc;
      glowLE.NeIRI(k,:) = ionoLE.Ne;
      glowLE.A5577(k,:) = ionoLE.A5577;
      glowLE.A4278(k,:) = ionoLE.A4278;
      glowLE.R5577(k,:) = ionoLE.R5577;
      glowLE.R4278(k,:) = ionoLE.R4278;
      glowLE.A6300(k,:) = ionoLE.A6300;
      glowLE.R6300(k,:) = ionoLE.R6300;
      
      ionoHE = calculate_rayleigh(time(iTime), glat, glon, Ap, energyBin, numFlux(iTime,:), energyBinGlow, filterHE);
      glowHE.Ne(k,:) = ionoHE.NeCalc;
      glowHE.NeIRI(k,:) = ionoHE.Ne;
      glowHE.A5577(k,:) = ionoHE.A5577;
      glowHE.A4278(k,:) = ionoHE.A4278;
      glowHE.R5577(k,:) = ionoHE.R5577;
      glowHE.R4278(k,:) = ionoHE.R4278;
      glowHE.A6300(k,:) = ionoHE.A6300;
      glowHE.R6300(k,:) = ionoHE.R6300;
      
      ionoSE = calculate_rayleigh(time(iTime), glat, glon, Ap, energyBin, numFlux(iTime,:), energyBinGlow, filterSE);
      glowSE.Ne(k,:) = ionoSE.NeCalc;
      glowSE.NeIRI(k,:) = ionoSE.Ne;
      glowSE.A5577(k,:) = ionoSE.A5577;
      glowSE.A4278(k,:) = ionoSE.A4278;
      glowSE.R5577(k,:) = ionoSE.R5577;
      glowSE.R4278(k,:) = ionoSE.R4278;
      glowSE.A6300(k,:) = ionoSE.A6300;
      glowSE.R6300(k,:) = ionoSE.R6300;
      
      iono = calculate_rayleigh(time(iTime), glat, glon, Ap, energyBin, zeros(1,length(energyBin)), energyBinGlow, filter);
     
      glow0.Ne(k,:) = iono.NeCalc;
      glow0.NeIRI(k,:) = iono.Ne;
      glow0.A5577(k,:) = iono.A5577;
      glow0.A4278(k,:) = iono.A4278;
      glow0.R5577(k,:) = iono.R5577;
      glow0.R4278(k,:) = iono.R4278;
      glow0.A6300(k,:) = iono.A6300;
      glow0.R6300(k,:) = iono.R6300;
      
      k = k+1;
end
%%
k = 1;
multiWaitbar('Progress Zero Energy Flux',0);
timeIndx = find_time(time,timeMinStr):1:find_time(time,timeMaxStr);
nTime = length(timeIndx);
dt = 1./nTime;
for iTime = timeIndx
    multiWaitbar('Progress Zero Energy Flux','Increment',dt);

    
    k = k+1;
end
glow0.time = time(timeIndx);
%%
glowSE.time = time(timeIndx);
glowHE.time = time(timeIndx);
glowLE.time = time(timeIndx);
glow.time = time(timeIndx);
toc

%%
X = dascData.latitude;
Y = dascData.longitude;
el = dascData.elevation;
az = dascData.azimuth;
Indx = ~isnan(X);
[I,J] = ndgrid(1:length(dascData.latitude),1:length(dascData.longitude));
Fi = scatteredInterpolant(X(Indx),Y(Indx),I(Indx),'nearest');
Fj = scatteredInterpolant(X(Indx),Y(Indx),J(Indx),'nearest');
vi = Fi(glat,glon);
vj = Fj(glat,glon);

[i,j] = get_neighbouring_pixels(vi,vj,az,el,1);

%%
i = squeeze(i);
j = squeeze(j);

[opticalSignal, timeOptical] = get_undersampled_asi(dascImage,i,j,dascData.time,timeMinStr,timeMaxStr);

%%

%%
figure;
plot(glowHE.A4278,iono.alt);
ylabel('km');
xlabel('Volume emission rates cm-3 s-1');
title('30-300 keV e^- precipitation');
%%
figure; 
plot(datetime(datevec(datestr(glowHE.time))),glowHE.A4278(:,26)); %90 km
hold on
plot(datetime(datevec(datestr(glowHE.time))),glowHE.A4278(:,107)); % 300 km
ylabel('4278A Volume emission rates cm-3 s-1');
title('30-300 keV e^- precipitation');
legend('90 km', '300 km');
%%
figure; 
plot(datetime(datevec(datestr(glow.time))),glowHE.R4278);
hold on;
plot(datetime(datevec(datestr(glow0.time))),glow0.R4278);
hold on;
title('>30 keV electron precipitation');
yyaxis right
ylabel('Counts');
opticalSignalFiltered = opticalSignal;
opticalSignalFiltered(opticalSignalFiltered>500)=nan;
plot(datetime(datevec(datestr(timeOptical))),opticalSignalFiltered,'Color',[0.5 0.5 0.5]);
legend('Glow - 4278A in R','No precipitation','DASC Counts','Location','NorthWest');
yyaxis left
ylabel('Rayleigh');
% ylim([550,1000]);

%%
figure; 

plot(datetime(datevec(datestr(glow.time))),glow.R5577-glow0.R5577,'g');
hold on;
plot(datetime(datevec(datestr(glow.time))),glow.R4278-glow0.R4278,'b');
hold on;
plot(datetime(datevec(datestr(glow.time))),glow.R6300-glow0.R6300,'r');
% hold on;
% plot(datetime(datevec(datestr(glow.time))),glow.R5577+glow.R4278+glow.R6300,'k');

title('Glow comparison with DASC data');
yyaxis right
ylabel('Counts');
opticalSignalFiltered = opticalSignal;
opticalSignalFiltered(opticalSignalFiltered>500)=nan;
plot(datetime(datevec(datestr(timeOptical))),opticalSignalFiltered,'Color',[0.5 0.5 0.5]);
legend('GLOW - 5577A','GLOW - 4278A','GLOW - 6300A','DASC Counts','Location','NorthWest');
yyaxis left
ylabel('Rayleigh');
% ylim([1000,12000]);

%%
figure; 
plot(datetime(datevec(datestr(glow.time))),glowHE.R5577-glow0.R5577,'g');
hold on;
plot(datetime(datevec(datestr(glow.time))),glowHE.R4278-glow0.R4278,'b');
hold on;
title('>30 keV electron precipitation');
% yyaxis right
% ylabel('Counts');
% opticalSignalFiltered = opticalSignal;
% opticalSignalFiltered(opticalSignalFiltered>500)=nan;
% plot(datetime(datevec(datestr(timeOptical))),opticalSignalFiltered,'Color',[0.5 0.5 0.5]);
legend('Glow - 5577A in R without airglow','Glow - 4278A in R without airglow','Location','NorthWest');
yyaxis left
ylabel('Rayleigh');

%%
totalPanelNo=2;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);

q1=p(1);
% 4a. Keogram, White light image
q1(1).select();
plot(glow.time,glowSE.R5577-glow0.R5577,'m');
hold on;
plot(glow.time,glowHE.R5577-glow0.R5577,'b');
ylim([0 500]);
[TTick,TTickLim]=label_time_axis(glow.time,false,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission','5577 Aº','[Rayleigh]'});
legend('>100 keV e^-','>30 keV e^-','Location','northwest');
% set(gca,'YScale','log');
% ylim([10,4000]);

q1(2).select();
plot(glow.time,glowSE.R4278-glow0.R4278,'m');
hold on;
plot(glow.time,glowHE.R4278-glow0.R4278,'b');
ylim([0 500]);

[TTick,TTickLim]=label_time_axis(glow.time,true,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission','4278 Aº','[Rayleigh]'});
legend('>100 keV e^-','>30 keV e^-','Location','northwest');


%%
totalPanelNo=3;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);

q=p(1);
% 4a. Keogram, White light image
q(2).select();
plot(glow.time,glowHE.R6300-glow0.R6300,'r');
hold on;
plot(glow.time,glowHE.R5577-glow0.R5577,'g');
hold on;
plot(glow.time,glowHE.R4278-glow0.R4278,'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,false,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','>30 keV e^-','[Rayleigh]'});
legend('6300Aº','5577Aº','4278Aº','Location','northwest');
ylim([0 400]);

q(3).select();
plot(glow.time,glowSE.R6300-glow0.R6300,'r');
hold on;
plot(glow.time,glowSE.R5577-glow0.R5577,'g');
hold on;
plot(glow.time,glowSE.R4278-glow0.R4278,'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,true,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','>100 keV e^-','[Rayleigh]'});
legend('6300Aº','5577Aº','4278Aº','Location','northwest');
ylim([0 100]);

q(1).select();
plot(glow.time,glow.R6300-glow0.R6300,'r');
hold on;
plot(glow.time,glow.R5577-glow0.R5577,'g');
hold on;
plot(glow.time,glow.R4278-glow0.R4278,'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,false,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','>0.05 keV e^-','[Rayleigh]'});
legend('6300Aº','5577Aº','4278Aº','Location','northwest');
ylim([0 7000]);

%%
totalPanelNo=3;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);

q=p(1);
% 4a. Keogram, White light image
q(2).select();
plot(glow.time,100.*(glowHE.R6300-glow0.R6300)./(glow.R6300-glow0.R6300),'r');
hold on;
plot(glow.time,100.*(glowHE.R5577-glow0.R5577)./(glow.R5577-glow0.R5577),'g');
hold on;
plot(glow.time,100.*(glowHE.R4278-glow0.R4278)./(glow.R4278-glow0.R4278),'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,false,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','>30 keV e^-','[%]'});
ylim([0 100]);
legend('6300Aº','5577Aº','4278Aº','Location','northwest');

q(3).select();
plot(glow.time,100.*(glowSE.R6300-glow0.R6300)./(glow.R6300-glow0.R6300),'r');
hold on;
plot(glow.time,100.*(glowSE.R5577-glow0.R5577)./(glow.R5577-glow0.R5577),'g');
hold on;
plot(glow.time,100.*(glowSE.R4278-glow0.R4278)./(glow.R4278-glow0.R4278),'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,true,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','>100 keV e^-','[%]'});
legend('6300Aº','5577Aº','4278Aº','Location','northwest');

q(1).select();
plot(glow.time,100.*(glowLE.R6300-glow0.R6300)./(glow.R6300-glow0.R6300),'r');
hold on;
plot(glow.time,100.*(glowLE.R5577-glow0.R5577)./(glow.R5577-glow0.R5577),'g');
hold on;
plot(glow.time,100.*(glowLE.R4278-glow0.R4278)./(glow.R4278-glow0.R4278),'b');
hold on;
[TTick,TTickLim]=label_time_axis(glow.time,false,1/6,timeMinStr,timeMaxStr);
ylabel({'Emission from','<30 keV e^-','[%]'});
legend('6300Aº','5577Aº','4278Aº','Location','northwest');
%%
figure; 
plot(datetime(datevec(datestr(glowLE.time))),glowLE.R4278); 
hold on; 
plot(datetime(datevec(datestr(glowHE.time))),glowHE.R4278);


%%
figure; 
loglog(ionoLE.energyBin,ionoLE.phitop,'r-'); 
hold on; 
loglog(ionoHE.energyBin,ionoHE.phitop,'b-'); 
hold on;
loglog(energyBin, numFlux(iTime,:),'k--');
hold on;

%%
figure;
semilogx(glow.Ne,iono.alt,'r');
hold on;
semilogx(electronDensity(:,iTime),pfisrNeAlt,'k');
ylim([60,200]);

%%
figure;
semilogx(glow.A5577,ionoLE.alt,'g');
hold on;
semilogx(glow.A4278,ionoLE.alt,'b');
%%


function iono = calculate_rayleigh(time, glat, glon, Ap, energyBin, energyFlux, energyBinGlow, filter)
    
    f=150;
    numFluxGlow = interp1(energyBin,energyFlux,energyBinGlow);
    numFluxGlow(isnan(numFluxGlow)) = 0;
    numFluxGlow = numFluxGlow.*filter; 
    iono = glowenergy(time, glat, glon, f, f, f, Ap, energyBinGlow, numFluxGlow);
    
end

function ener=egrid(eMax, nBins)
   
    ener = zeros(1,nBins);
    for n=1:21
       ener(n) = 0.5.*n;
    end
    ener(n+1:end) = logspace(log10(11),log10(eMax),nBins-21);

    del(1) = 0.5;

    for n=2:nBins
     del(n) = ener(n)-ener(n-1);
    end
  
    for n=1:nBins
     ener(n) = ener(n) - del(n)/2.0;
    end
  
end
function [i_neighbourhood,j_neighbourhood] = get_neighbouring_pixels(vi,vj,asiaz,asiel,beamWidth)
[I, J] = ndgrid(1:length(asiaz),1:length(asiel));
i_neighbourhood = zeros(size(vi,1),size(vi,2),64);
j_neighbourhood = zeros(size(vi,1),size(vi,2),64);
asizn = 90-asiel;
for i = 1:1:size(vi,1)
    for j= 1:1:size(vi,2)
        az0 = asiaz(vi(i,j),vj(i,j));
        zn0 = asizn(vi(i,j),vj(i,j));
        % translation of a circle in polar coordinates
        indx = asizn.^2 + zn0.^2 - 2.*asizn.*zn0.*cosd(asiaz-az0)<=beamWidth.*0.5;
        i_neighbourhood(i,j,1:length(I(indx))) = I(indx);
        j_neighbourhood(i,j,1:length(I(indx))) = J(indx);
    end
end

end

function [asi_undersampled,time_undersampled] = get_undersampled_asi(asi,I,J,time,timeMinStr,timeMaxStr)
% asi  is the all sky image [nTime x nPixels x nPixels]
% I are indices of interest [nBeams x nAlt or nEnergy x col number]
% J are indices of interest [nBeams x nAlt or nEnergy x row number]
% I(p,q,:) - the corresponding pixels in asi are to be averaged

timeIndx = find_time(time,timeMinStr):1:find_time(time,timeMaxStr);

% Get linear index numbers for a slice of ASI in time
nSize = size(squeeze(asi(1,:,:)));
nTime = size(asi,1);
asi_undersampled = zeros(length(timeIndx),1);

linearInd.arr = calculate_linear_index(nSize, I, J); 

% Slicing ASI in time
multiWaitbar('UndersampleASI',0);
dt = 1./length(timeIndx);
k = 1;
for iT = timeIndx
    multiWaitbar('UndersampleASI','Increment',dt);
    ASI = squeeze(asi(iT,:,:));
    asi_undersampled(k) = mean(ASI(linearInd.arr));    
    k = k+1;
end
time_undersampled = time(timeIndx);
end

function linearInd = calculate_linear_index (matrixSize, i_n, j_n)

linearInd = squeeze(sub2ind(matrixSize,i_n(i_n(:)>0),j_n(j_n(:)>0)));

end