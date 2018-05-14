function [dataASILastDay] = create_DASC_hdf5(localStoreDir,outputH5FileStr,...
    projectionAltitude,minElevation,...
    minTimeStr,maxTimeStr,...
    calFileAz,calFileEl, setDownloadDASCFlag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<9
    setDownloadDASCFlag = true;
end
if nargin<8 || isempty(calFileAz)
    calFileAz = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '26Mar2008',filesep,'PKR_Cal_before_2011',filesep,'PKR_20111006_AZ_10deg.FITS'];
end
if nargin<7 || isempty(calFileEl)
    calFileEl = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '26Mar2008',filesep,'PKR_Cal_before_2011',filesep,'PKR_20111006_EL_10deg.FITS'];
end
if nargin<6
    maxTimeStr=[];
end
if nargin<5
    minElevation=30;
end
if nargin<4
    projectionAltitude=110;
end
if nargin<3
    minTimeStr=[];
end
if nargin<2
    outputH5FileStr='outputDASC.h5';
end
if nargin<1 
    localStoreDir = [intialize_root_path,'LargeFiles',filesep,'DASC',filesep];
end

imageSize=512;

if ~isempty(maxTimeStr) && ~isempty(minTimeStr)
    ndays = days(datenum(datestr(maxTimeStr,'dd-mmm-yyyy'))- ...
        datenum(datestr(minTimeStr,'dd-mmm-yyyy')));
    dayArray = datenum(datestr(minTimeStr,'dd-mmm-yyyy')):...
        1:datenum(datestr(minTimeStr,'dd-mmm-yyyy'))+ndays;
elseif ~isempty(maxTimeStr) && isempty(minTimeStr)
    dayArray = datenum(dateStr(maxTimeStr,'dd-mmm-yyyy'));
elseif ~isempty(minTimeStr) && isempty(maxTimeStr)
    dayArray = datenum(dateStr(minTimeStr,'dd-mmm-yyyy'));
else
    error('No date is specified to process DASC data');
end
multiWaitbar('Creating DASC HDF5 File',0);
for idays=1:1:length(dayArray)
    if setDownloadDASCFlag == true
        fprintf(['Status: Downloading FITS files ',datestr(dayArray(idays),'dd-mmm-yyyy')]);
        download_DASC_FITS(datestr(dayArray(idays)),localStoreDir);
    end
    
    localDASCDirPath(idays,:) = ([localStoreDir,filesep,datestr(dayArray(idays),'yyyymmdd')]); 
    
    fileStr = get_files_in_folder(localDASCDirPath(idays,:));
    fitsTimeStamp = fitsfiletimestamp(fileStr);
    timeASI = (fitsTimeStamp-datenum('jan-01-1970'))*(24*3600);
    timeASI = unix_to_matlab_time(timeASI);
    
    nTimeASI = length(timeASI);
    timeStartIndx = 1; timeEndIndx = nTimeASI;
    DASCCalFile.az = calFileAz; DASCCalFile.el = calFileEl; 
    azOldRes = fitsread(DASCCalFile.az);
    elOldRes = fitsread(DASCCalFile.el);
    data.ASI = nan(nTimeASI,imageSize.^2);
    data.lat = nan(nTimeASI,imageSize.^2);
    data.lon = nan(nTimeASI,imageSize.^2);
    data.az_new = nan(nTimeASI,imageSize.^2);
    data.el_new = nan(nTimeASI,imageSize.^2);
    data.timeDASC = nan(nTimeASI,1);
%     sensorloc = nan(3,1);
    multiWaitbar('Reading data',0);
    for itime=timeStartIndx:1:timeEndIndx
         multiWaitbar('Reading data','Increment',1/timeEndIndx);
        ASIDataStr = strcat(localDASCDirPath(idays,:),filesep,(fileStr(itime)));
        try
            [ASI(1,:), lat(1,:), lon(1,:),...
                az_new(1,:), el_new(1,:),...
                sensorloc, timeDASC] =...
                DASC_aer_to_geodetic(char(ASIDataStr),azOldRes,elOldRes,...
                imageSize,minElevation,projectionAltitude);
            notNANCoords = 1:length(lat);
            data.ASI(itime,notNANCoords) = ASI; 
            data.lat(itime,notNANCoords) = lat;
            data.lon(itime,notNANCoords) = lon;
            data.az_new(itime,notNANCoords) = az_new;
            data.el_new(itime,notNANCoords) = el_new;
            data.timeDASC(itime,1) = timeDASC;
            message(itime) = strcat('Successfully loaded: ',ASIDataStr);
        catch ME
            message(itime)=strcat('Could not load file: ',ASIDataStr);
            data.timeDASC(itime,1) = timeASI(itime);
        end
        
    end
    
    %% Writing Data
    fieldNames = fieldnames(data);
    nFields = length(fieldNames);
    multiWaitbar('Writing data',0);

%     if ~isfile(outputH5FileStr)
%         fid=HF5.create(outputH5FileStr);
%         H5F.close(fid);
%     end
        
    for iField = 1:1:nFields
    multiWaitbar('Writing data','Increment',1/nFields);
        if iField <6
    h5create(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        size(data.(fieldNames{iField})'),'ChunkSize',[50 80],'Deflate',9);
    h5write(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        (data.(fieldNames{iField}))');
    h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        'Dimensions','nTime x nCoords');
        else
    h5create(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        size(data.(fieldNames{iField})'),'ChunkSize',[1 50],'Deflate',9);
    h5write(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        (data.(fieldNames{iField}))');
    h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',fieldNames{iField}],...
        'Dimensions','nTime x 1');
        end
            
    end
h5create(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','azCalData'],...
    size(azOldRes'),'ChunkSize',[50 80],'Deflate',9);
h5write(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','azCalData'],...
    (azOldRes'));
h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','azCalData'],...
    'Calibration File Used',calFileAz);

h5create(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','elCalData'],...
    size(elOldRes'),'ChunkSize',[50 80],'Deflate',9);
h5write(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','elCalData'],...
    (elOldRes'));
h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','elCalData'],...
    'Calibration File Used',calFileEl);

dset = (message)';
dset_details.Location = ['/DASC/',datestr(dayArray(idays),'yyyymmdd')];
dset_details.Name = 'message';

attr = 'nTime x nMessageLength';
attr_details.Name = 'Dimensions';
attr_details.AttachedTo = ['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/','message'];
attr_details.AttachType = 'dataset';

hdf5write(outputH5FileStr,dset_details,dset,...
    attr_details,attr,'WriteMode','append');

dset = sensorloc;
dset_details.Name = 'sensorloc';
attr = '3x1 [lat,lon,alt-meters]';
attr_details.AttachedTo = ['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/',dset_details.Name];
attr_details.AttachType = 'dataset';
hdf5write(outputH5FileStr,dset_details,dset,...
    attr_details,attr,'WriteMode','append');

h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/'],...
    'creation_date',datestr(now));
h5writeatt(outputH5FileStr,['/DASC/',datestr(dayArray(idays),'yyyymmdd'),'/'],...
    'duration_contained_in_file',[datestr(data.timeDASC(1)),' - ',datestr(data.timeDASC(end))]);

multiWaitbar('Creating DASC HDF5 File','Increment',1/length(dayArray));    
multiWaitbar('Writing data','Reset');
multiWaitbar('Reading data','Reset');

end

multiWaitbar('Creating DASC HDF5 File','Reset');
dataASILastDay = data;
end