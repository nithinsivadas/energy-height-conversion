%% Energy spectra

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/PFISR_002_006/Nithin/Data/';
    storeDir = '/media/nithin/PFISR_002_006/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
else
%     error(['Not configured for this computer: ',get_computer_name]);
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
end

%% Getting input files
fileNameList = struct2cell(dir([storeDir,'*_pfisrData.h5']));
filePathStr = strcat(storeDir,string(fileNameList(1,:)'));
%% Calculating and writing conductivity back into input files
energyBin = logspace(3,6,25)';
% for i=1:1:length(filePathStr)
for i=1:1:length(filePathStr)
    fileName = filePathStr(i);
    [data] = calculate_energy(fileName,energyBin);
    write_energy(data,fileName);
end


function [data] = calculate_energy(fileName,energyBin)
   
   stepsHour = hours(1); %Calculates conductivity and writes it in steps of x hours
   row = @(x,y) strcmp(string(x.Path),y);
   in = read_h5_data(fileName);
   alt = in.Data{row(in,'/alt')}';
   electronDensity = in.Data{row(in,'/inputData/Ne')};
   dNeFrac = in.Data{row(in,'/inputData/dNeFrac')}; % dNeFrac is actually dNe
   latitude = in.Data{row(in,'/site/latitude')};
   longitude = in.Data{row(in,'/site/longitude')};
   time = in.Data{row(in,'/time')};
   ttime = datetime(time,'ConvertFrom','datenum');

    
%    mode = 0;
%    outputs = [];
%    setPlotConductivity = false;
   
   hourArr = ttime(1):stepsHour:ttime(end);
   nhour = length(hourArr);

   for iHour = 1:1:nhour
       minTimeIndx = find_time(time,datenum(hourArr(iHour)));
       if iHour ~=nhour 
        maxTimeIndx = find_time(time,datenum(hourArr(iHour+1)));
       else
        maxTimeIndx = length(time);
       end
       medTimeIndx = round((minTimeIndx + maxTimeIndx)/2);
       lat1=latitude*ones(size(alt));
       lon1=longitude*ones(size(alt));
       A = get_energy_dep_matrix(alt,energyBin,lat1,lon1,time(medTimeIndx));
       
       timeIndx = minTimeIndx:1:maxTimeIndx;
       
       if length(timeIndx)==1
       data.energyFlux(timeIndx,:) = nan;
        data.denergyFlux(timeIndx,:) = nan;
        data.MSE(timeIndx,:) = nan;
        data.qInverted(timeIndx,:) = nan;
        data.q(timeIndx,:) = nan;
        data.dq(timeIndx,:) = nan;
       else
       [dq, q] = get_error_in_q(electronDensity(timeIndx,:)',...
            dNeFrac(timeIndx,:)',alt,time(timeIndx)',2);
        dataInv = get_inverted_flux(q,dq,time(timeIndx),alt,energyBin,A);
        data.energyFlux(timeIndx,:) = dataInv.energyFlux';
        data.denergyFlux(timeIndx,:) = dataInv.dEnergyFlux';
        data.MSE(timeIndx,:) = dataInv.MSE';
        data.qInverted(timeIndx,:) = dataInv.qInverted;
        data.q(timeIndx,:) = q;
        data.dq(timeIndx,:) = dq;       
       end
        
        
        data.time(timeIndx) = time(timeIndx);
       
   end
   data.energyBin = energyBin;
   data.altitude = alt;
end

function write_energy(data,fileName)
    write_h5_dataset(fileName,'/energy/energyFlux',data.energyFlux,1);
    write_h5_dataset(fileName,'/energy/dEnergyFlux',data.denergyFlux,1);
    write_h5_dataset(fileName,'/energy/energyBin',data.energyBin,0);
    write_h5_dataset(fileName,'/energy/MSE',data.MSE,1);
    write_h5_dataset(fileName,'/energy/q',data.q,1);
    write_h5_dataset(fileName,'/energy/dq',data.dq,1);
    write_h5_dataset(fileName,'/energy/qInverted',data.qInverted,1);
end

