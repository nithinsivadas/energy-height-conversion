%% Generating themis spacecraft data
clear all;
multiWaitbar('CloseAll');
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

% Output Path
rootPath = 'G:\My Drive\Research\Projects\Yiqun Yu\Data\ASCII\Processed\';
%% Collecting data into a structure
tic
probes={'tha','thd','the'};
nProbes = length(probes);
iP = 1./nProbes;
multiWaitbar('1.Probes...',0);
for iProbe = 1:1:nProbes
    % Main Time (is SST time)
    time = dataPAD.(probes{iProbe}).sst.time(2:end);
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
    multiWaitbar('2.Calculating Loss-cone',0);
    dt = 1./nTime;
    for iTime = 1:1:nTime
        [padData.(probes{iProbe}).lossConeAngle(iTime,1)] =...
            find_loss_cone_angle(temp.(probes{iProbe}).XYZ_GEO(1),...
            temp.(probes{iProbe}).XYZ_GEO(2),...
            temp.(probes{iProbe}).XYZ_GEO(3),...
            time(iTime),magFieldNo,thisMaginput(iTime,:),300,85);
        multiWaitbar('2.Calculating Loss-cone','Increment',dt);
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
    padData.(probes{iProbe}).Units.lcDiffEflux = 'eV cm-2 sr-1 s-1 eV-1';
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

    multiWaitbar('3.Calculating Loss-cone-Flux',0);
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
            
            laSr=(convert_deg_to_sr(thisLcArray(end)));
            
            padData.(probes{iProbe}).lcEfluxYi(iTime,iEnergy) = lcEfluxYi;
            padData.(probes{iProbe}).lcDiffEfluxYi(iTime,iEnergy) = lcEfluxYi./laSr;
            padData.(probes{iProbe}).fitYi.x0(iTime,iEnergy,:) = yYi;
            padData.(probes{iProbe}).fitYi.MSE(iTime,iEnergy) = MSEYi;
            
            padData.(probes{iProbe}).lcEfluxAr(iTime,iEnergy) = lcEfluxAr;
            padData.(probes{iProbe}).lcDiffEfluxAr(iTime,iEnergy) = lcEfluxAr./laSr;
            padData.(probes{iProbe}).fitAr.x0(iTime,iEnergy,:) = yAr;
            padData.(probes{iProbe}).fitAr.MSE(iTime,iEnergy) = MSEAr;
            
            padData.(probes{iProbe}).lcEfluxLi(iTime,iEnergy) = lcEfluxLi;
            padData.(probes{iProbe}).lcDiffEfluxLi(iTime,iEnergy) = lcEfluxLi./laSr;
        end
    multiWaitbar('3.Calculating Loss-cone-Flux','Increment',dt);    
    end

    %% Writing File
    multiWaitbar('4.Writing file',0);
    lossConeFluxEstimatorStr = 'Li';
    fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
    write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);
    multiWaitbar('4.Writing file',0.3);
    lossConeFluxEstimatorStr = 'Ar';
    fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
    write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);
    multiWaitbar('4.Writing file',0.6);
    lossConeFluxEstimatorStr = 'Yi';
    fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr);
    write_data(fileID,padData,probes,iProbe,lossConeFluxEstimatorStr);
    multiWaitbar('4.Writing file',1);
    
    %% Setting MultiWaitbar to 0, and incrementing
    multiWaitbar('1.Probes...','Increment',iP);
    multiWaitbar('2.Calculating Loss-cone',0);
    multiWaitbar('3.Calculating Loss-cone-Flux',0);
    multiWaitbar('4.Writing file',0);
end
toc

function write_data(fileID, padData, probes, iProbe, lossConeFluxEstimatorStr)
    
   if strcmp(lossConeFluxEstimatorStr,'Ar')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionAr);
        lcEflux = padData.(probes{iProbe}).lcDiffEfluxAr;
    elseif strcmp(lossConeFluxEstimatorStr,'Yi')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionYi);
        lcEflux = padData.(probes{iProbe}).lcDiffEfluxYi;
    elseif strcmp(lossConeFluxEstimatorStr,'Li')
        funcStr = 'Linear Extrapolation';
        lcEflux = padData.(probes{iProbe}).lcDiffEfluxLi;
    else
        error('Unknown lossConeFluxEstimatorStr');
    end
    
    % Correcting for loss-cone-flux which have negative values
    lcEflux(lcEflux<0) = 0;
    % Correcting for Bx,By,Bz-flux which have negative values
    padData.(probes{iProbe}).Bx_GSE(abs(padData.(probes{iProbe}).Bx_GSE)>10^5)=nan;
    padData.(probes{iProbe}).By_GSE(abs(padData.(probes{iProbe}).By_GSE)>10^5)=nan;
    padData.(probes{iProbe}).Bz_GSE(abs(padData.(probes{iProbe}).Bz_GSE)>10^5)=nan;
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
        '# DataFormat: [%s] \n'...
        '# \n'],...
        probes{iProbe},...
        padData.(probes{iProbe}).Description.lossConeAngle,...
        funcStr,...
        dataFormat);
    
    % Header
    fprintf(fileID,headerFormat,'DateTime','UT','X_GSE','Y_GSE','Z_GSE','Bx_GSE',...
        'By_GSE','Bz_GSE','alpha_LC',padData.(probes{iProbe}).energyBin');
    
    % Units
    fprintf(fileID,unitHeaderFormat,'dd-MMM-yyyy HH:mm:ss','Posix[s]','[RE]','[RE]','[RE]','[nT]',...
        '[nT]','[nT]','[deg]',string(repmat('eV/cm2 sr s eV',42,1)));
   
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
    
    fileName =  ['lc_',probes{iProbe},'_',lossConeFluxEstimatorStr,'_diffEflux_nithin.dat'];
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
    try
        [y, resnormy] = lsqcurvefit(padFun,x0,pa,padData,[],[],options);
        MSE = resnormy/(length(pa)-2);
        padFunNew = @(x) padFun(y,x);
        lossConeEnergyFlux = trapz(convert_deg_to_sr(lcArray),padFunNew(lcArray)); %[eV/cm2 s eV]
    catch ME;
        y = [nan nan];
        MSE = nan;
        lossConeEnergyFlux = nan;
    end
end
function [lossConeEnergyFlux] = get_loss_cone_flux_linear_interp(pa,padData,lcArray)
    try
        padNew = interp1(pa,padData,lcArray,'linear','extrap');
        lossConeEnergyFlux = trapz(convert_deg_to_sr(lcArray),padNew);
    catch ME;
        lossConeEnergyFlux = nan;
    end
    
end
function sr = convert_deg_to_sr(angle)
    sr = 4*pi*(sind(angle./2)).^2;
end