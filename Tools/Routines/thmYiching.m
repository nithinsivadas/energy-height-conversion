%% Generating themis spacecraft data

% Load data
if ispc
    load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\thd_data_26_Mar_2008.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\pitchAngleDistribution20080326.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\magneticField20080326.mat')
%     load ('C:\Users\nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\thd_state.mat');
else
    load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
%     load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/thd_state.mat');
end
[maginput,timeMaginput]=generate_maginput('C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5','26 Mar 2008 8:00','26 Mar 2008 13:00');
[stateData, matFilePath] = process_themis_data('26-Mar-2008', [initialize_root_path,'LargeFiles',filesep],...
    'thd','state');

%% Time
timeStr = '26-Mar-2008 11:10';
energy = 50; %in KeV
%% Pitch angle distribution
% paFluxSST = thd.pitchAngleDistribution.sst.flux; %Differential energy flux
% time = thd.pitchAngleDistribution.sst.time;
% pa = cell2mat(thd.pitchAngleDistribution.sst.pitchAngle);
paFluxSST = dataPAD.thd.sst.paFlux;
paHighRes = 0:1:180;
energyBinIndx = find_altitude(dataPAD.thd.sst.energyBin,energy*1000);
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
time = dataPAD.thd.sst.time;
thisTimeIndx = find_time(time,timeStr);

pad = squeeze(paFluxSST(thisTimeIndx,energyBinIndx,:));
pa = dataPAD.thd.sst.pa;
figure;
plot(pa,pad);

% pad = @(x,xdata) x(1)*(sind(xdata)).^x(2);
padAmertyev = @(x,xdata) x(1)*(cosd(xdata)).^2 + x(2)*(sind(xdata)).^2;

% x0 = [10^5,10];
x0Amertyev = [10^3 10^3];

% [x, resnorm] = lsqcurvefit(pad,x0,pa,paFluxSST(thisTimeIndx,:));
[y, resnormy] = lsqcurvefit(padAmertyev,x0Amertyev,pa,pad);
% hold on;
% plot(paHighRes,pad(x,paHighRes),'-r');
hold on;
plot(paHighRes,padAmertyev(y,paHighRes),'-r');
set(gca,'YScale','log');
legend('Data','Fit');

% chisquared=resnorm/(length(pa)-2)
disp(['Chi-square: ',num2str(resnormy/(length(pa)-2))]);
disp(['Energy: ',num2str(dataPAD.thd.sst.energyBin(energyBinIndx)/1000),' keV']);

%% Check if PAD from toshi's files adds up to enregy flux to confirm units
energyFluxSST = thd.energySpectra.energyFlux(1,32:end); %SST Flux
paEnergyFluxSST = squeeze(paFluxSST(2,:,:));
paSr = 4*pi*(sind(pa./2)).^2;
integratedEnergyFlux = (trapz(paSr,paEnergyFluxSST')/(4*pi))*10^4; % solid angle maximum is 4pi Sr
% Therefore paEnergyFlux is eV/cm^2 s sr eV
figure;plot(dataPAD.thd.sst.energyBin,integratedEnergyFlux); hold on;
plot(dataPAD.thd.sst.energyBin,energyFluxSST);
set(gca,'XScale','log');
set(gca,'YScale','log');
legend('Integrated Flux over PA','Energy Flux from SPEDAS');
%% Process themis data
stateTimeIndx = find_time(stateData.thd.state.time,timeStr);
xGEO = stateData.thd.state.XYZ_GEO(stateTimeIndx,:);
xGSE = stateData.thd.state.XYZ_GSE(stateTimeIndx,:);
thisMaginput = interp1(timeMaginput,maginput,stateData.thd.state.time(stateTimeIndx));
[la,alt] = find_loss_cone_angle(xGEO(1),xGEO(2),xGEO(3),...
    stateData.thd.state.time(stateTimeIndx),4,thisMaginput,300,85);

%% Calculating the loss-cone flux
% Normalizing the pitch angle distribution


function [la,altitude]=find_loss_cone_angle(x1,x2,x3,thisTime,magneticFieldModelNo,maginput,nPoints,stopAlt)
alphaArray=linspace(0,3,nPoints);
    for i=1:1:length(alphaArray)
        [~,~,xGEO] = onera_desp_lib_find_mirror_point(magneticFieldModelNo,...
        [0,0,0,0,0],1,thisTime,x1,x2,x3,alphaArray(i),maginput);
        R(i) = sqrt(xGEO(1).^2+xGEO(2).^2+xGEO(3).^2);
    end
C = define_universal_constants;
altitude=((R-1)*C.RE)/1000;
la = interp1(altitude(~isnan(altitude)),alphaArray(~isnan(altitude)),stopAlt);
end