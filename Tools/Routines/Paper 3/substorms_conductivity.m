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

fileName = filePathStr(1);
[data,input] = calculate_conductivity(fileName);



function [data,input] = calculate_conductivity(fileName)
    tic
   in = read_h5_data(fileName);
   alt = in.Data{1}';
   electronDensity = in.Data{8}';
   latitude = in.Data{11};
   longitude = in.Data{12};
   time = in.Data{7}';
%    minAlt = 80;
%    maxAlt = max(alt);
%    electronDensity = crop_altitude(electronDensity,alt,minAlt,maxAlt);
%    alt = crop_altitude(alt,alt,minAlt,maxAlt);
   mode = 1;
   outputs = 'all';
   setPlotConductivity = false;
   iriData = iri2016f90(time(1), alt, latitude, longitude);
   alt1=alt(:)';
   lat1=latitude*ones(size(alt1));
   lon1=longitude*ones(size(alt1));
   time1=time(1)*ones(size(alt1));
   msisData = msis_irbem(time1, [alt1',lat1',lon1']);
   toc
   tic
   [data, input] = get_conductivity_v2( alt, electronDensity(:,1),...
    latitude, longitude, time(1), mode, outputs,setPlotConductivity,...
    iriData, msisData);   
   toc
    
end

