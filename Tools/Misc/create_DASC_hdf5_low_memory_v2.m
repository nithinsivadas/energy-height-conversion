function [dataASILastDay] = create_DASC_hdf5_low_memory_v2(localStoreDir,outputH5FileStr,...
    projectionAltitudeDASC,minElevation,...
    minTimeStr,maxTimeStr,...
    calFileLat,calFileLon, setDownloadDASCFlag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<9
    setDownloadDASCFlag = true;
end
if nargin<8 || isempty(calFileLon)
    calFileLon = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glon.dat'];
if nargin<7 || isempty(calFileLat)
    calFileLat = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glat.dat'];
end

end
if nargin<6
    maxTimeStr=[];
end
if nargin<5
    minTimeStr=[];
end
if nargin<4
    minElevation=30;
end
if nargin<3 || isempty(projectionAltitudeDASC)
    projectionAltitudeDASC=110;
end
if nargin<2
    outputH5FileStr='outputDASC.h5';
end
if nargin<1 
    localStoreDir = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep];
end

imageSize=1024;
fprintf(['Projection Alt of DASC ',num2str(projectionAltitudeDASC),' km']);

if ~isempty(maxTimeStr) && ~isempty(minTimeStr)
    ndays = days(datenum(datestr(maxTimeStr,'dd-mmm-yyyy'))- ...
        datenum(datestr(minTimeStr,'dd-mmm-yyyy')));
    dayArray = datenum(datestr(minTimeStr,'dd-mmm-yyyy')):...
        1:datenum(datestr(minTimeStr,'dd-mmm-yyyy'))+ndays;
elseif ~isempty(maxTimeStr) && isempty(minTimeStr)
    dayArray = datenum(datestr(maxTimeStr,'dd-mmm-yyyy'));
elseif ~isempty(minTimeStr) && isempty(maxTimeStr)
    dayArray = datenum(datestr(minTimeStr,'dd-mmm-yyyy'));
else
    error('No date is specified to process DASC data');
end
multiWaitbar('Creating DASC HDF5 File',0);
start = [1 1 1];
start5 = [1 1];
saveMemoryFlag = isunix;

fileLatID = fopen(calFileLat);
latCal = cell2mat(textscan(fileLatID,repmat('%8.4f',1,1024)));
latCal=interp_nans(latCal)';
fclose(fileLatID);

fileLonID = fopen(calFileLon);
lonCal = cell2mat(textscan(fileLonID,repmat('%8.4f',1,1024)));
lonCal=interp_nans(lonCal)';
fclose(fileLonID);

sensorloc=[65.1260,-147.4789,689 ]; %poker flat location

[azCal,elCal]=geodetic2aer(latCal,lonCal,110,sensorloc(1),sensorloc(2),sensorloc(3)/1000,wgs84Ellipsoid('km'));

for idays=1:1:length(dayArray)
    if setDownloadDASCFlag == true
        fprintf(['Status: Downloading FITS files ',datestr(dayArray(idays),'dd-mmm-yyyy')]);
        download_DASC_FITS(datestr(dayArray(idays)),localStoreDir);
    end
    
    localDASCDirPath(idays,:) = ([localStoreDir,datestr(dayArray(idays),'yyyymmdd')]); 
    
    fileStr = get_files_in_folder(localDASCDirPath(idays,:),'*.FITS');
    timeASI = fitsfiletimestamp(fileStr);
    
    if idays==1 && ~isempty(minTimeStr)
        [fileStr, timeASI]=crop_time(reshape(fileStr,[],1),timeASI,datenum(minTimeStr), timeASI(end));
    end
    if idays==length(dayArray) && ~isempty(maxTimeStr)
        [fileStr, timeASI]=crop_time(reshape(fileStr,[],1),timeASI,timeASI(1),datenum(maxTimeStr));
    end
    
    nTimeASITotal = length(timeASI);
    dkTime = 800; % To allow for low memory
    nkTime = ceil(nTimeASITotal/dkTime);
    multiWaitbar('Reading data',0);
    multiWaitbar('Writing data',0);
    data.lat = nan(imageSize,imageSize);
    data.lon = nan(imageSize,imageSize);
    data.alt = nan(imageSize,imageSize);
    data.az = nan(imageSize,imageSize);
    data.el = nan(imageSize,imageSize);
    data.lowAzGradientFilter = nan(imageSize,imageSize);
    data.minElFilter = nan(imageSize,imageSize);
    sizeData1 = [imageSize,imageSize];
    
    for kTime=1:nkTime    
        timeEndIndx = min(kTime*dkTime,nTimeASITotal);
        timeStartIndx = 1 + (kTime-1)*dkTime;
        nTimeASI = length(timeStartIndx:timeEndIndx);
%         timeStartIndx = 1 ; timeEndIndx = nTimeASI;
        
        ASI = nan(nTimeASI,imageSize,imageSize); %size = [nTimeASI, imageSize.^2];
        data.time = nan(nTimeASI,1); %size = [nTimeASI, 1]

        fieldNames = fieldnames(data);
        nFields = length(fieldNames);
        sizeData = [imageSize,imageSize, nTimeASI];
        sizeData5 = [1, nTimeASI];

        count = [imageSize imageSize nTimeASI];
        count5 = [1 nTimeASI];
        if(sizeData(1)<64)chunkSize(1,1)=sizeData(1);else chunkSize(1,1)=64;end
        if(sizeData(2)<64)chunkSize(1,2)=sizeData(2);else chunkSize(1,2)=64;end
        if(sizeData(3)<80)chunkSize(1,3)=sizeData(3);else chunkSize(1,3)=80;end
        if(sizeData1(1)<50)chunkSize1(1,1)=sizeData1(1);else chunkSize1(1,1)=50;end
        if(sizeData1(2)<80)chunkSize1(1,2)=sizeData1(2);else chunkSize1(1,2)=80;end  
        if(sizeData5(1)<50)chunkSize5(1,1)=sizeData5(1);else chunkSize5(1,1)=50;end
        if(sizeData5(2)<80)chunkSize5(1,2)=sizeData5(2);else chunkSize5(1,2)=80;end  

        if idays==1 && kTime==1
            h5create(outputH5FileStr,['/DASC/','ASI'],...
                [sizeData(1) sizeData(2) Inf],'ChunkSize',chunkSize,'Deflate',9);
            h5writeatt(outputH5FileStr,['/DASC/','ASI'],...
                'Dimensions','nTime x nImageSize x nImageSize');
            for iField = 1:1:nFields
                if iField < nFields 
                    h5create(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        sizeData1,'ChunkSize',chunkSize1,'Deflate',9);
                    h5writeatt(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        'Dimensions','nImageSize x nImageSize');
                else

                    h5create(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        [sizeData5(1) Inf],'ChunkSize',chunkSize5,'Deflate',9);
                    h5writeatt(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        'Dimensions','nTime x 1');
                    h5writeatt(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        'Units','Posix time');
                end
            end
         data.az = azCal;
         data.el = elCal;
         data.lat = latCal;
         data.lon = lonCal;
         data.alt = projectionAltitudeDASC.*ones(size(data.el));
         data.minElFilter = double(elCal<minElevation);
        end
        
        itempTime=1;
        for itime=timeStartIndx:1:timeEndIndx

            multiWaitbar('Reading data','Increment',1/nTimeASITotal);
            ASIDataStr = strcat(localDASCDirPath(idays,:),filesep,(fileStr(itime)));

            try
                [ASI(itempTime,:,:),data.time(itempTime,1)]=read_DASC_fits(char(ASIDataStr),imageSize);
                message(itempTime) = strcat('Successfully loaded: ',ASIDataStr);
            catch ME
                message(itempTime)=strcat('Could not load file: ',ASIDataStr,' because of: ', ME.message, 'in Function ', ME.stack.file);
                data.time(itempTime,1) = timeASI(itime);
            end
        itempTime = itempTime+1;
        end

        %% Writing Data

        data.time = posixtime(datetime(data.time,'ConvertFrom','datenum'));
        if idays==1 && kTime == 1
            for iField = 1:1:nFields-1
                h5write(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        (data.(fieldNames{iField}))');  
            end
        multiWaitbar('Writing data','Increment',0.2/nkTime);
        h5create(outputH5FileStr,['/DASC/','azCalData'],...
            size(azCal'),'ChunkSize',[50 80],'Deflate',9);
        h5write(outputH5FileStr,['/DASC/','azCalData'],...
            (azCal'));
        h5writeatt(outputH5FileStr,['/DASC/','azCalData'],...
            'Calibration File Used',calFileLat);

        h5create(outputH5FileStr,['/DASC/','elCalData'],...
            size(elCal'),'ChunkSize',[50 80],'Deflate',9);
        h5write(outputH5FileStr,['/DASC/','elCalData'],...
            (elCal'));
        h5writeatt(outputH5FileStr,['/DASC/','elCalData'],...
            'Calibration File Used',calFileLon);
        end
        multiWaitbar('Writing data','Increment',0.4/nkTime);
        %Writing 'time'
        iField = nFields; 
        h5write(outputH5FileStr,['/DASC/',fieldNames{iField}],...
                        (data.(fieldNames{iField}))',start5,count5);
        multiWaitbar('Writing data','Increment',0.5/nkTime);
        %Writing image
        h5write(outputH5FileStr,['/DASC/','ASI'],...
                    permute(ASI,[3 2 1]),start,count);       
        multiWaitbar('Writing data','Value',kTime*1./nkTime);
        % Writing error messages
        messageTotal(start5(2):1:start5(2)+nTimeASI-1) = message;
        multiWaitbar('Creating DASC HDF5 File','Increment',1/(length(dayArray)*kTime));    
        clearvars message;
        start = start + [0 0 nTimeASI];
        start5 = start5 + [0 nTimeASI];
    end
     if idays==1
        startTime = timeASI(1);
     end
     if idays==length(dayArray)
         endTime = timeASI(end);
     end
     multiWaitbar('Reading data','Reset');
     multiWaitbar('Writing data','Reset');
end

dset = (messageTotal)';
dset_details.Location = '/DASC/';
dset_details.Name = 'message';

attr = 'nTime x nMessageLength';
attr_details.Name = 'Dimensions';
attr_details.AttachedTo = ['/DASC/','message'];
attr_details.AttachType = 'dataset';

hdf5write(outputH5FileStr,dset_details,dset,...
    attr_details,attr,'WriteMode','append');

dset = sensorloc;
dset_details.Name = 'sensorloc';
attr = '3x1 [lat,lon,alt-meters]';
attr_details.AttachedTo = ['/DASC/',dset_details.Name];
attr_details.AttachType = 'dataset';
hdf5write(outputH5FileStr,dset_details,dset,...
    attr_details,attr,'WriteMode','append');

h5writeatt(outputH5FileStr,'/DASC/',...
    'creation_date',datestr(now));
h5writeatt(outputH5FileStr,'/DASC/',...
    'duration_contained_in_file',[datestr(startTime),' - ',datestr(endTime)]);

add_ASI_background_to_hdf5('dasc',outputH5FileStr); % Adding the backgroud

multiWaitbar('Creating DASC HDF5 File','Reset');
dataASILastDay = data;
end