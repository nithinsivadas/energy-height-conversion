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
for i=1:1:length(filePathStr)
    fileName = filePathStr(i);
    [data] = calculate_conductivity(fileName);
    write_conductivity(data,fileName);
end


function [data] = calculate_conductivity(fileName)
   
   stepsHour = hours(1); %Calculates conductivity and writes it in steps of x hours
   row = @(x,y) strcmp(string(x.Path),y);
   in = read_h5_data(fileName);
   alt = in.Data{row(in,'/alt')}';
   electronDensity = in.Data{row(in,'/inputData/Ne')};
   latitude = in.Data{row(in,'/site/latitude')};
   longitude = in.Data{row(in,'/site/longitude')};
   time = in.Data{row(in,'/time')};
   ttime = datetime(time,'ConvertFrom','datenum');

    
   mode = 0;
   outputs = [];
   setPlotConductivity = false;
   
   hourArr = ttime(1):stepsHour:ttime(end);
   nhour = length(hourArr);
   alt1=alt(:)';
   lat1=latitude*ones(size(alt1));
   lon1=longitude*ones(size(alt1));
   for iHour = 1:1:nhour
       minTimeIndx = find_time(time,datenum(hourArr(iHour)));
       if iHour ~=nhour 
        maxTimeIndx = find_time(time,datenum(hourArr(iHour+1)));
       else
        maxTimeIndx = length(time);
       end
       medTimeIndx = round((minTimeIndx + maxTimeIndx)/2);
       iriData = iri2016f90(time(medTimeIndx), alt, latitude, longitude);
       time1=time(medTimeIndx)*ones(size(alt1));
       msisData = msis_irbem(time1, [alt1',lat1',lon1']);
       for iTime = minTimeIndx:1:maxTimeIndx
        [tempData] = get_conductivity_v2( alt, electronDensity(iTime,:)',...
        latitude, longitude, time(iTime), mode, outputs,setPlotConductivity,...
        iriData, msisData);
        data.pedersonConductivity(iTime,:) = tempData.pedersonConductivity;
        data.hallConductivity(iTime,:) = tempData.hallConductivity;
        data.time(iTime) = time(iTime);
       end
   end
   data.altitude = tempData.altitude;
end

function write_conductivity(data,fileName)
    write_h5_dataset(fileName,'/conductivity/hall',real(data.hallConductivity),1);
    write_h5_dataset(fileName,'/conductivity/pederson',real(data.pedersonConductivity),1);
    write_h5_dataset(fileName,'/conductivity/ihall',imag(data.hallConductivity),1);
    write_h5_dataset(fileName,'/conductivity/ipederson',imag(data.pedersonConductivity),1);
end

