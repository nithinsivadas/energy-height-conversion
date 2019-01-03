%% Generating themis spacecraft data
%% Need to add a function to find the trapped particles!
clear all;
multiWaitbar('CloseAll');
%% Load data
if ispc
%     load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\thd_data_26_Mar_2008.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\pitchAngleDistribution20080326.mat')
%     load('G:\My Drive\Research\Projects\Yiqun Yu\Data\magneticField20080326.mat')
    load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
%     load ('C:\Users\nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\thd_state.mat');
else
    load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
%     load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/thd_state.mat');
end

% Output Path
rootPath = 'G:\My Drive\Research\Projects\Yiqun Yu\Data\ASCII\Processed\';
[stateData, matFilePath] = process_themis_data('26-Mar-2008', [initialize_root_path,'LargeFiles',filesep],...
    'tha,thd,the','state');
%% Collecting data into a structure
tic
probes={'tha','thd','the'};
nProbes = length(probes);
iP = 1./nProbes;
multiWaitbar('1.Probes...',0);
for iProbe = 1:1:nProbes
        
    %% Calculating Trapped Flux
    time = padData.(probes{iProbe}).time;
    nTime = length(padData.(probes{iProbe}).time);
    nEnergyBin=length(padData.(probes{iProbe}).energyBin);
    %pitch angle distrubtion
    nPitchAngle = length(padData.(probes{iProbe}).pitchAngle);
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
    
    padFunArtemyev = @(x,xdata) x(:,1).*(cosd(xdata)).^2 + x(:,2).*(sind(xdata)).^2; %AV Artemyev - ?2014
    padFunYi = @(x,xdata) x(1)*(sind(xdata)).^x(2);
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
    pa = padData.(probes{iProbe}).pitchAngle;
    padData.(probes{iProbe}).pitchAngleDistributionFunctionYi = padFunYi;
    padData.(probes{iProbe}).pitchAngleDistributionFunctionAr = padFunArtemyev;
    padData.(probes{iProbe}).Units.trappedEflux = 'eV/cm2 s eV';
    padData.(probes{iProbe}).Units.trappedDiffEflux = 'eV/cm2 sr s eV';
    padData.(probes{iProbe}).Units.trappedNflux = '#/cm2 s eV';
    padData.(probes{iProbe}).Units.trappedDiffNflux = '#/cm2 s sr eV';
    padData.(probes{iProbe}).Units.pa = 'deg';
    padData.(probes{iProbe}).Description.pitchAngle = 'Pitch Angle Bins';
    padData.(probes{iProbe}).Description.energyBin = 'Energy bins from ESA & SST instruments combined';
    padData.(probes{iProbe}).Description.X_GSE = 'X_GSE coordinate of THEMIS probe in RE';
    padData.(probes{iProbe}).Description.Y_GSE = 'Y_GSE coordinate of THEMIS probe in RE';
    padData.(probes{iProbe}).Description.Z_GSE = 'Z_GSE coordinate of THEMIS probe in RE';
    padData.(probes{iProbe}).Description.Bx_GSE = 'Bx_GSE local spin-averaged magnetic field in nT';
    padData.(probes{iProbe}).Description.By_GSE = 'By_GSE local spin-averaged magnetic field in nT';
    padData.(probes{iProbe}).Description.Bz_GSE = 'Bz_GSE local spin-averaged magnetic field in nT';

    multiWaitbar('3.Calculating Trapped Flux',0);
    dt = 1./nTime;
    for iTime = 1:1:nTime
        thisTrappedArray = linspace(80,110,50);
        for iEnergy = 1:1:nEnergyBin
            thisPad = squeeze(temp.(probes{iProbe}).pad(iTime,iEnergy,:));
            med = median(thisPad,1);
            x0Yi = [med 1];
            x0Ar = [med med];
            [trappedEfluxYi,yYi,MSEYi] = get_trapped_flux(padFunYi,x0Yi,pa,thisPad,thisTrappedArray,options);
            [trappedEfluxAr,yAr,MSEAr] = get_trapped_flux(padFunArtemyev,x0Ar,pa,thisPad,thisTrappedArray,options);
            [trappedEfluxLi] = get_trapped_flux_linear_interp(pa,thisPad,thisTrappedArray);
            
            trappedSr=(convert_deg_to_sr(thisTrappedArray(end))-convert_deg_to_sr(thisTrappedArray(1)));
            
            padData.(probes{iProbe}).trappedEfluxYi(iTime,iEnergy) = trappedEfluxYi;
            padData.(probes{iProbe}).trappedDiffEfluxYi(iTime,iEnergy) = trappedEfluxYi./trappedSr;
            padData.(probes{iProbe}).fitYi.x0(iTime,iEnergy,:) = yYi;
            padData.(probes{iProbe}).fitYi.MSE(iTime,iEnergy) = MSEYi;
            
            padData.(probes{iProbe}).trappedEfluxAr(iTime,iEnergy) = trappedEfluxAr;
            padData.(probes{iProbe}).trappedDiffEfluxAr(iTime,iEnergy) = trappedEfluxAr./trappedSr;
            padData.(probes{iProbe}).fitAr.x0(iTime,iEnergy,:) = yAr;
            padData.(probes{iProbe}).fitAr.MSE(iTime,iEnergy) = MSEAr;
            
            padData.(probes{iProbe}).trappedEfluxLi(iTime,iEnergy) = trappedEfluxLi;
            padData.(probes{iProbe}).trappedDiffEfluxLi(iTime,iEnergy) = trappedEfluxLi./trappedSr;
        end
    padData.(probes{iProbe}).trappedNfluxYi(iTime,:) = padData.(probes{iProbe}).trappedEfluxYi(iTime,:)./padData.(probes{iProbe}).energyBin';
    padData.(probes{iProbe}).trappedDiffNfluxYi(iTime,:) = padData.(probes{iProbe}).trappedDiffEfluxYi(iTime,:)./padData.(probes{iProbe}).energyBin';
    padData.(probes{iProbe}).trappedNfluxAr(iTime,:) = padData.(probes{iProbe}).trappedEfluxAr(iTime,:)./padData.(probes{iProbe}).energyBin';
    padData.(probes{iProbe}).trappedDiffNfluxAr(iTime,:) = padData.(probes{iProbe}).trappedDiffEfluxAr(iTime,:)./padData.(probes{iProbe}).energyBin';
    padData.(probes{iProbe}).trappedNfluxLi(iTime,:) = padData.(probes{iProbe}).trappedEfluxLi(iTime,:)./padData.(probes{iProbe}).energyBin';
    padData.(probes{iProbe}).trappedDiffNfluxLi(iTime,:) = padData.(probes{iProbe}).trappedDiffEfluxLi(iTime,:)./padData.(probes{iProbe}).energyBin';
    multiWaitbar('3.Calculating Trapped-Flux','Increment',dt);    
    end

      
    %% Setting MultiWaitbar to 0, and incrementing
    multiWaitbar('1.Probes...','Increment',iP);
    multiWaitbar('2.Calculating Loss-cone',0);
    multiWaitbar('3.Calculating Loss-cone-Flux',0);
end

toc

function write_data(fileID, padData, probes, iProbe, lossConeFluxEstimatorStr,lossConeFluxStr)
    
   if strcmp(lossConeFluxEstimatorStr,'Ar')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionAr);
    elseif strcmp(lossConeFluxEstimatorStr,'Yi')
        funcStr = func2str(padData.(probes{iProbe}).pitchAngleDistributionFunctionYi);
    elseif strcmp(lossConeFluxEstimatorStr,'Li')
        funcStr = 'Linear Extrapolation';
    else
        error('Unknown lossConeFluxEstimatorStr');
   end

    lcflux = padData.(probes{iProbe}).([lossConeFluxStr,lossConeFluxEstimatorStr]);
    
    % Correcting for loss-cone-flux which have negative values
    lcflux(lcflux<0) = 0;
    % Correcting for Bx,By,Bz-flux which have negative values
    padData.(probes{iProbe}).Bx_GSE(abs(padData.(probes{iProbe}).Bx_GSE)>10^5)=nan;
    padData.(probes{iProbe}).By_GSE(abs(padData.(probes{iProbe}).By_GSE)>10^5)=nan;
    padData.(probes{iProbe}).Bz_GSE(abs(padData.(probes{iProbe}).Bz_GSE)>10^5)=nan;
    % Defining data output format
    headerFormat = ['%20s ',repmat('%10s ',1,8+5),repmat('%15.3f ',1,42),'\n'];
    unitHeaderFormat = ['%20s ',repmat('%10s ',1,8+5),repmat('%15s ',1,42),'\n'];
    dataFormat = ['%20s %10.0f ',repmat('%10.7f ',1,3),repmat('%10.3f ',1,3+5),...
        '%10.4f ',repmat('%15.8e ',1,42),'\n'];
    
    % Comment
    fprintf(fileID,['# Loss cone energy flux estimated from THEMIS ESA & SST Measurements \n',...
        '# Probe: %s\n',...
        '# %s\n',...
        '# Input parameter of magnetic field model taken from omni database\n',...
        '# Loss cone flux estimated using fit function: %s\n'...
        '# Loss cone flux values < 0 have been set to 0\n'...
        '# %s\n'...
        '# DataFormat: [%s] \n'...
        '# \n'],...
        probes{iProbe},...
        padData.(probes{iProbe}).Description.lossConeAngle,...
        funcStr,...
        padData.(probes{iProbe}).Description.magFieldStrFP,...
        dataFormat);
    
    % Header
    fprintf(fileID,headerFormat,'DateTime','UT','X_GSE','Y_GSE','Z_GSE',...
        'Lat_FP','Lon_FP','MLat_FP','MLon_FP','MLT_FP',...
        'Bx_GSE','By_GSE','Bz_GSE','alpha_LC',padData.(probes{iProbe}).energyBin');
    
    % Units
    fprintf(fileID,unitHeaderFormat,'dd-MMM-yyyy HH:mm:ss','Posix[s]','[RE]','[RE]','[RE]',...
        'deg','deg','deg','deg','hr',...
        '[nT]','[nT]','[nT]','[deg]',string(repmat(padData.(probes{iProbe}).Units.lcDiffNflux,42,1)));
   
    % Data 
    fprintf(fileID,dataFormat,[string(padData.(probes{iProbe}).timeStr),...
        padData.(probes{iProbe}).UTtime,...
        padData.(probes{iProbe}).XYZ_GSE(:,1),...
        padData.(probes{iProbe}).XYZ_GSE(:,2),...
        padData.(probes{iProbe}).XYZ_GSE(:,3),...
        padData.(probes{iProbe}).latFoot,...
        padData.(probes{iProbe}).lonFoot,...
        padData.(probes{iProbe}).mlatFoot,...
        padData.(probes{iProbe}).mlonFoot,...
        padData.(probes{iProbe}).mltFoot,...
        padData.(probes{iProbe}).Bx_GSE,...
        padData.(probes{iProbe}).By_GSE,...
        padData.(probes{iProbe}).Bz_GSE,...
        padData.(probes{iProbe}).lossConeAngle,...
        lcflux]');
    
    %Close file
    fclose(fileID);
end

function fileID = create_file(rootPath,probes,iProbe,lossConeFluxEstimatorStr,lossConeFluxStr)
    
%     fileName =  ['lc_',probes{iProbe},'_',lossConeFluxEstimatorStr,'_Nflux_nithin.dat'];
    fileName =  ['lc_',probes{iProbe},'_',lossConeFluxEstimatorStr,'_',lossConeFluxStr(3:end),'_with_foot_point.dat'];
    filePath = [rootPath,fileName];
    fileID = fopen(filePath,'w');
    
end

function [la,altitude]=find_loss_cone_angle(x1,x2,x3,thisTime,magneticFieldModelNo,maginput,nPoints,stopAlt)
    alphaArray=linspace(0,10,nPoints);
    [~,~,xGEO] = onera_desp_lib_find_mirror_point(magneticFieldModelNo,...
                 [0,0,0,0,0],1,thisTime,x1,x2,x3,alphaArray,maginput);
    R = sqrt(xGEO(:,1).^2+xGEO(:,2).^2+xGEO(:,3).^2);
    
    C = define_universal_constants;
    altitude=((R-1)*C.RE)/1000;
    try
        la = interp1(altitude(~isnan(altitude)),alphaArray(~isnan(altitude)),stopAlt,'linear','extrap');
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

function [trappedEnergyFlux,y,MSE] = get_trapped_flux(padFun,x0,pa,padData,trapedArray,options)
    try
        [y, resnormy] = lsqcurvefit(padFun,x0,pa,padData,[],[],options);
        MSE = resnormy/(length(pa)-2);
        padFunNew = @(x) padFun(y,x);
        trappedEnergyFlux = trapz(convert_deg_to_sr(trapedArray),padFunNew(trapedArray)); %[eV/cm2 s eV]
    catch ME;
        y = [nan nan];
        MSE = nan;
        trappedEnergyFlux = nan;
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

function [trappedEnergyFlux] = get_trapped_flux_linear_interp(pa,padData,trappedArray)
    try
        padNew = interp1(pa,padData,trappedArray,'linear','extrap');
        trappedEnergyFlux = trapz(convert_deg_to_sr(trappedArray),padNew);
    catch ME;
        trappedEnergyFlux = nan;
    end
end

function sr = convert_deg_to_sr(angle)
    sr = 4*pi*(sind(angle./2)).^2;
end