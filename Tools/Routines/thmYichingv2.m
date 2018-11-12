%% Generating themis spacecraft data
clear all;

%% Load data
if ispc
%     load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\thd_data_26_Mar_2008.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\pitchAngleDistribution20080326.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\magneticField20080326.mat')
%     load ('C:\Users\nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\thd_state.mat');
else
    load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
%     load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/thd_state.mat');
end
[maginput,timeMaginput]=generate_maginput('C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5','26 Mar 2008 8:00','26 Mar 2008 13:00');
[stateData, matFilePath] = process_themis_data('26-Mar-2008', [initialize_root_path,'LargeFiles',filesep],...
    'tha,thd,the','state');

%% Collecting data into a structure
probes={'tha','thd','the'};

iProbe = 2;

% Main Time (is SST time)
time = dataPAD.(probes{iProbe}).sst.time(2:5);
padData.(probes{iProbe}).time = time;
padData.(probes{iProbe}).timeStr = datestr(time,'dd-mmm-yyyy HH:MM:SS');
padData.(probes{iProbe}).UTtime = posixtime(datetime(padData.(probes{iProbe}).timeStr,'InputFormat','dd-MMM-yyyy HH:mm:ss'));
nTime = length(time);
%energyBin
padData.(probes{iProbe}).energyBin = [dataPAD.(probes{iProbe}).esa.energyBin;...
    dataPAD.(probes{iProbe}).sst.energyBin];
nEnergyBin=length(padData.(probes{iProbe}).energyBin);
%pitch angle
padData.(probes{iProbe}).pitchAngle = dataPAD.(probes{iProbe}).sst.pa;
nPitchAngle = length(padData.(probes{iProbe}).pitchAngle);
%pitch angle distrubtion
temp.(probes{iProbe}).pad = zeros(nTime,nEnergyBin,nPitchAngle);
temp.(probes{iProbe}).pad(1:1:nTime,1:1:31,1:1:nPitchAngle) =...
    interp1(dataPAD.(probes{iProbe}).esa.time,dataPAD.(probes{iProbe}).esa.paFlux,...
    time,'linear');
temp.(probes{iProbe}).pad(1:1:nTime,32:nEnergyBin,1:1:nPitchAngle) =...
    interp1(dataPAD.(probes{iProbe}).sst.time,dataPAD.(probes{iProbe}).sst.paFlux,...
        time,'linear');    
%state of probe
temp.(probes{iProbe}).XYZ_GEO =...
    interp1(stateData.(probes{iProbe}).state.time,stateData.(probes{iProbe}).state.XYZ_GEO,...
    time,'linear');
padData.(probes{iProbe}).XYZ_GSE =...
    interp1(stateData.(probes{iProbe}).state.time,stateData.(probes{iProbe}).state.XYZ_GSE,...
    time,'linear');
%magnetic field
padData.(probes{iProbe}).Bx_GSE =...
    interp1(dataMag.(probes{iProbe}).time,dataMag.(probes{iProbe}).BxGSE,...
    time,'linear');
padData.(probes{iProbe}).By_GSE =...
    interp1(dataMag.(probes{iProbe}).time,dataMag.(probes{iProbe}).ByGSE,...
    time,'linear');
padData.(probes{iProbe}).Bz_GSE =...
    interp1(dataMag.(probes{iProbe}).time,dataMag.(probes{iProbe}).BzGSE,...
    time,'linear');
%% Loss cone angle 
% loss cone angle is considered not varying with energy (since the
% variation is approx 0.05 deg from 1keV to 100 keV

thisMaginput = interp1(timeMaginput',maginput,time);
magFieldNo = 4;
multiWaitbar('Calculating Loss-cone',0);
dt = 1./nTime;
for iTime = 1:1:nTime
    multiWaitbar('Calculating Loss-cone','Increment',dt);
    [padData.(probes{iProbe}).lossConeAngle(iTime,1)] =...
        find_loss_cone_angle(temp.(probes{iProbe}).XYZ_GEO(1),...
        temp.(probes{iProbe}).XYZ_GEO(2),...
        temp.(probes{iProbe}).XYZ_GEO(3),...
        time(iTime),magFieldNo,thisMaginput(iTime,:),300,85);
end

%% Calculating Loss-cone Flux
padFunArtemyev = @(x,xdata) x(:,1).*(cosd(xdata)).^2 + x(:,2).*(sind(xdata)).^2; %AV Artemyev - ?2014
padFunYi = @(x,xdata) x(1)*(sind(xdata)).^x(2);
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
pa = padData.(probes{iProbe}).pitchAngle;
la = padData.(probes{iProbe}).lossConeAngle;
padData.(probes{iProbe}).pitchAngleDistributionFunctionYi = padFunYi;
padData.(probes{iProbe}).pitchAngleDistributionFunctionAr = padFunArtemyev;
padData.(probes{iProbe}).Units.lcEflux = 'eV cm-2 s-1 eV-1';
padData.(probes{iProbe}).Units.la = 'deg';
padData.(probes{iProbe}).Units.pa = 'deg';
padData.(probes{iProbe}).Description.lossConeAngle = ['Loss cone angle, calulated using: ',...
    find_irbem_magFieldModelStr(magFieldNo),' Model'];
padData.(probes{iProbe}).Description.pa = 'Pitch Angle Bins';
padData.(probes{iProbe}).Description.energyBin = 'Energy bins from ESA & SST instruments combined';
padData.(probes{iProbe}).Description.X_GSE = 'X_GSE coordinate of THEMIS probe in RE';
padData.(probes{iProbe}).Description.Y_GSE = 'Y_GSE coordinate of THEMIS probe in RE';
padData.(probes{iProbe}).Description.Z_GSE = 'Z_GSE coordinate of THEMIS probe in RE';
padData.(probes{iProbe}).Description.Bx_GSE = 'Bx_GSE local spin-averaged magnetic field in nT';
padData.(probes{iProbe}).Description.By_GSE = 'By_GSE local spin-averaged magnetic field in nT';
padData.(probes{iProbe}).Description.Bz_GSE = 'Bz_GSE local spin-averaged magnetic field in nT';

for iTime = 1:1:nTime
    thisLcArray = linspace(0.0001,la(iTime),50);
    for iEnergy = 1:1:nEnergyBin
        thisPad = squeeze(temp.(probes{iProbe}).pad(iTime,iEnergy,:));
        med = median(thisPad,1);
        x0Yi = [med 1];
        x0Ar = [med med];
        [lcEfluxYi,yYi,MSEYi] = get_loss_cone_flux(padFunYi,x0Yi,pa,thisPad,thisLcArray,options);
        [lcEfluxAr,yAr,MSEAr] = get_loss_cone_flux(padFunArtemyev,x0Ar,pa,thisPad,thisLcArray,options);
        [lcEfluxLi] = get_loss_cone_flux_linear_interp(pa,thisPad,thisLcArray);
        padData.(probes{iProbe}).lcEfluxYi(iTime,iEnergy) = lcEfluxYi;
        padData.(probes{iProbe}).fitYi.x0(iTime,iEnergy,:) = yYi;
        padData.(probes{iProbe}).fitYi.MSE(iTime,iEnergy) = MSEYi;
        padData.(probes{iProbe}).lcEfluxAr(iTime,iEnergy) = lcEfluxAr;
        padData.(probes{iProbe}).fitAr.x0(iTime,iEnergy,:) = yAr;
        padData.(probes{iProbe}).fitAr.MSE(iTime,iEnergy) = MSEAr;
        padData.(probes{iProbe}).lcEfluxLi(iTime,iEnergy) = lcEfluxLi;
        
    end
end



%% Sample Printing
rootPath = 'G:\My Drive\Research\Projects\Yiqun Yu\Data\ASCII\Processed\';
lossConeFluxEstimatorStr = 'Li';
fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);
lossConeFluxEstimatorStr = 'Ar';
fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);
lossConeFluxEstimatorStr = 'Yi';
fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);

%% Displaying sample
% fprintf('\n');
% write_data([],padData,probes,iProbe,'Ar');
% fprintf('\n');
% write_data([],padData,probes,iProbe,'Yi');

% %% Sample
% iTime = 2;
% iEnergy = 30;
% thisPad = squeeze(temp.(probes{iProbe}).pad(iTime,iEnergy,:));
% thisLcArray = linspace(0.0001,la(iTime),50);
% med = median(thisPad,1);
% x0Yi = [med 1];
% x0Ar = [med med];
% [lcEfluxYi,yYi,MSEYi] = get_loss_cone_flux(padFunYi,x0Yi,pa,thisPad,thisLcArray,options);
% [lcEfluxAr,yAr,MSEAr] = get_loss_cone_flux(padFunArtemyev,x0Ar,pa,thisPad,thisLcArray,options);
% [lcEfluxLi] = get_loss_cone_flux_linear_interp(pa,thisPad,thisLcArray);
% 
% padHighRes = linspace(0,180,300);
% figure;
% plot(pa,thisPad,'.-k'); hold on;
% plot(padHighRes,padFunArtemyev(yAr,(padHighRes)),'r');
% hold on;
% plot(padHighRes,padFunYi(yYi,(padHighRes)),'g');
% set(gca,'XScale','linear');
% set(gca,'YScale','log');
% legend('Data','Artemyev Fit','Yiching Fit');
% 
% fprintf(['\nloss cone energy flux (lc angle= %3.2f deg, Enegry = %5.2f keV): \n',...
%     ' Yi      : %15.8e eV/cm^2 s eV \n',...
%     ' Artemyev: %15.8e eV/cm^2 s eV \n',...
%     ' Linear  : %15.8e eV/cm^2 s eV \n'],la(iTime),padData.(probes{iProbe}).energyBin(iEnergy)/1000,lcEfluxYi,lcEfluxAr,lcEfluxLi);
% 
% 
% %% Time
% timeStr = '26-Mar-2008 11:10';
% energy = 50; %in KeV
% %% Pitch angle distribution
% % paFluxSST = thd.pitchAngleDistribution.sst.flux; %Differential energy flux
% % time = thd.pitchAngleDistribution.sst.time;
% % pa = cell2mat(thd.pitchAngleDistribution.sst.pitchAngle);
% paFluxSST = dataPAD.thd.sst.paFlux;
% paHighRes = 0:1:180;
% energyBinIndx = find_altitude(dataPAD.thd.sst.energyBin,energy*1000);
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% time = dataPAD.thd.sst.time;
% thisTimeIndx = find_time(time,timeStr);
% 
% pad = squeeze(paFluxSST(thisTimeIndx,energyBinIndx,:));
% pa = dataPAD.thd.sst.pa;
% figure;
% plot(pa,pad);
% 
% % pad = @(x,xdata) x(1)*(sind(xdata)).^x(2);
% padAmertyev = @(x,xdata) x(1)*(cosd(xdata)).^2 + x(2)*(sind(xdata)).^2;
% 
% % x0 = [10^5,10];
% x0Amertyev = [10^3 10^3];
% 
% % [x, resnorm] = lsqcurvefit(pad,x0,pa,paFluxSST(thisTimeIndx,:));
% [y, resnormy] = lsqcurvefit(padAmertyev,x0Amertyev,pa,pad);
% % hold on;
% % plot(paHighRes,pad(x,paHighRes),'-r');
% hold on;
% plot(paHighRes,padAmertyev(y,paHighRes),'-r');
% set(gca,'YScale','log');
% legend('Data','Fit');
% 
% % chisquared=resnorm/(length(pa)-2)
% disp(['Chi-square: ',num2str(resnormy/(length(pa)-2))]);
% disp(['Energy: ',num2str(dataPAD.thd.sst.energyBin(energyBinIndx)/1000),' keV']);
% 
% %% Check if PAD from toshi's files adds up to enregy flux to confirm units
% energyFluxSST = thd.energySpectra.energyFlux(1,32:end); %SST Flux
% paEnergyFluxSST = squeeze(paFluxSST(2,:,:));
% paSr = 4*pi*(sind(pa./2)).^2;
% integratedEnergyFlux = (trapz(paSr,paEnergyFluxSST')/(4*pi))*10^4; % solid angle maximum is 4pi Sr
% % Therefore Toshi's file unit is [eV/cm^2 s sr eV]
% figure;plot(dataPAD.thd.sst.energyBin,integratedEnergyFlux); hold on;
% plot(dataPAD.thd.sst.energyBin,energyFluxSST);
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% legend('Integrated Flux over PA','Energy Flux from SPEDAS');
% %% Process themis data
% stateTimeIndx = find_time(stateData.thd.state.time,timeStr);
% xGEO = stateData.thd.state.XYZ_GEO(stateTimeIndx,:);
% xGSE = stateData.thd.state.XYZ_GSE(stateTimeIndx,:);
% thisMaginput = interp1(timeMaginput,maginput,stateData.thd.state.time(stateTimeIndx));
% [la,alt] = find_loss_cone_angle(xGEO(1),xGEO(2),xGEO(3),...
%     stateData.thd.state.time(stateTimeIndx),4,thisMaginput,300,85);
% 
% %% Calculating the loss-cone flux
% % Normalizing the pitch angle distribution
% 



function write_data(fileID, padData, probes, iProbe, lossConeFluxEstimatorStr)
    
   if strcmp(lossConeFluxEstimatorStr,'Ar')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionAr);
        lcEflux = padData.(probes{iProbe}).lcEfluxAr;
    elseif strcmp(lossConeFluxEstimatorStr,'Yi')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionYi);
        lcEflux = padData.(probes{iProbe}).lcEfluxYi;
    elseif strcmp(lossConeFluxEstimatorStr,'Li')
        funcStr = 'Linear Extrapolation';
        lcEflux = padData.(probes{iProbe}).lcEfluxLi;
    else
        error('Unknown lossConeFluxEstimatorStr');
    end
    
    % Correcting for loss-cone-flux which have negative values
    lcEflux(lcEflux<0) = 0;
    
    % Defining data output format
    headerFormat = ['%20s ',repmat('%10s ',1,8),repmat('%15.3f ',1,42),'\n'];
    unitHeaderFormat = ['%20s ',repmat('%10s ',1,8),repmat('%15s ',1,42),'\n'];
    dataFormat = ['%20s %10.0f ',repmat('%10.7f ',1,3),repmat('%10.3f ',1,3),...
        '%10.4f ',repmat('%15.8e ',1,42),'\n'];
    
    % Comment
    fprintf(fileID,['# Loss cone energy flux estimated from THEMIS ESA & SST Measurements \n',...
        '# Probe: %s\n',...
        '# %s\n',...
        '# Input parameter of magnetic field model taken from omni database\n',...
        '# Loss cone flux estiamted using fit function: %s\n'...
        '# Loss cone flux values < 0 have been set to 0\n'...
        '# \n'],...
        probes{iProbe},...
        padData.(probes{iProbe}).Description.la,...
        funcStr);
    
    % Header
    fprintf(fileID,headerFormat,'DateTime','UT','X_GSE','Y_GSE','Z_GSE','Bx_GSE',...
        'By_GSE','Bz_GSE','alpha_LC',padData.(probes{iProbe}).energyBin');
    
    % Units
    fprintf(fileID,unitHeaderFormat,'dd-MMM-yyyy HH:mm:ss','Posix[s]','[RE]','[RE]','[RE]','[nT]',...
        '[nT]','[nT]','[deg]',string(repmat('[eV]',42,1)));
   
    % Data 
    fprintf(fileID,dataFormat,[string(padData.(probes{iProbe}).timeStr),...
        padData.(probes{iProbe}).UTtime,...
        padData.(probes{iProbe}).XYZ_GSE(:,1),...
        padData.(probes{iProbe}).XYZ_GSE(:,2),...
        padData.(probes{iProbe}).XYZ_GSE(:,3),...
        padData.(probes{iProbe}).Bx_GSE,...
        padData.(probes{iProbe}).By_GSE,...
        padData.(probes{iProbe}).Bz_GSE,...
        padData.(probes{iProbe}).lossConeAngle,...
        lcEflux]');
    
    %Close file
    fclose(fileID);
end

function fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr)
    
    fileName =  ['lc_',probes{iProbe},'_',lossConeFluxEstimatorStr,'_nithin.dat'];
    filePath = [rootPath,fileName];
    fileID = fopen(filePath,'w');
    
end

function [la,altitude]=find_loss_cone_angle(x1,x2,x3,thisTime,magneticFieldModelNo,maginput,nPoints,stopAlt)
alphaArray=linspace(0,3,nPoints);
    for i=1:1:length(alphaArray)
        [~,~,xGEO] = onera_desp_lib_find_mirror_point(magneticFieldModelNo,...
        [0,0,0,0,0],1,thisTime,x1,x2,x3,alphaArray(i),maginput);
        R(i) = sqrt(xGEO(1).^2+xGEO(2).^2+xGEO(3).^2);
    end
C = define_universal_constants;
altitude=((R-1)*C.RE)/1000;
try
la = interp1(altitude(~isnan(altitude)),alphaArray(~isnan(altitude)),stopAlt);
catch ME;
    la = nan;
end
end

function [lossConeEnergyFlux,y,MSE] = get_loss_cone_flux(padFun,x0,pa,padData,lcArray,options)
    [y, resnormy] = lsqcurvefit(padFun,x0,pa,padData,[],[],options);
    MSE = resnormy/(length(pa)-2);
    padFunNew = @(x) padFun(y,x);
%     lax=linspace(0,la,50);
    lossConeEnergyFlux = trapz(convert_deg_to_sr(lcArray),padFunNew(lcArray)); %[eV/cm2 s eV]
%     lossConeEnergyFlux = integral(padFunNew,0,la); %eV/cm2 s eV
end
function [lossConeEnergyFlux] = get_loss_cone_flux_linear_interp(pa,padData,lcArray)
%     lax=linspace(0,la,50);
    padNew = interp1(pa,padData,lcArray,'linear','extrap');
    lossConeEnergyFlux = trapz(convert_deg_to_sr(lcArray),padNew);
%     lossConeEnergyFlux = integral(padFunNew,0,la); %eV/cm2 s eV
end
function sr = convert_deg_to_sr(angle)
    sr = 4*pi*(sind(angle./2)).^2;
end