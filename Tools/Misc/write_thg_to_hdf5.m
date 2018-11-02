function [input]=write_thg_to_hdf5(h5OutputFile,ASI,varargin)
    % write_thg_to_hdf5.m Writes thg data to hdf5 file under /sitename/. 
    % e.g. /GAKO/
    % The function will append data to the file if it already has data
    % under the group name. 

    p = inputParser;
    validationFcn = @(x) length(size(x))==3;
    
    addParameter(p,'lat',0);
    addParameter(p,'lon',0);
    addParameter(p,'az',0);
    addParameter(p,'el',0);
    addParameter(p,'alt',0);
    addParameter(p,'altIndx',2);
    addParameter(p,'mlat',0);
    addParameter(p,'mlon',0);
    addParameter(p,'time',0);
    addParameter(p,'sensorloc',0, @(x) prod(size(x))>=2 && prod(size(x))<4);
    addParameter(p,'siteCode','unknown');
    

    addRequired(p,'h5OutputFile',@(x)contains(x,{'.h5','.hdf5'})); 
    addRequired(p,'ASI',validationFcn);
   
    parse(p,h5OutputFile,ASI,varargin{:});
    
    input = p.Results;
    input.lat = squeeze(input.lat(:,:,input.altIndx));
    input.lon = squeeze(input.lon(:,:,input.altIndx));
    input.mlat = squeeze(input.mlat(:,:,input.altIndx));
    input.mlon = squeeze(input.mlon(:,:,input.altIndx));
    input.alt = repmat(input.alt(input.altIndx)/1000,size(input.az));
    fieldNames = fieldnames(input);
    nFields = length(fieldNames);
    
    imageSize = size(ASI,1);
    nTimeASI = length(p.Results.time);
    start = [1 1 1];
    start5 = [1 1];
    
    sizeData = [imageSize,imageSize, nTimeASI];
    sizeData1 = [imageSize,imageSize];
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
    
    groupName = ['/',upper(input.siteCode),'/'];
    
    % Creating dataset if previously non-existent
    [datasetExists, datasetInfo]=ish5dataset(h5OutputFile,groupName);
    if ~datasetExists
        h5create(h5OutputFile,[groupName,'ASI'],...
                [Inf Inf Inf],'ChunkSize',chunkSize,'Deflate',9);
        h5writeatt(h5OutputFile,[groupName,'ASI'],...
            'Dimensions','nTime x nImageSize x nImageSize');
        
        [~,fields2D] =  intersect(fieldNames,{'lat','lon','az','el',...
            'mlat','mlon','alt'});
        for iField = fields2D'
            h5create(h5OutputFile,[groupName,fieldNames{iField}],...
                sizeData1,'ChunkSize',chunkSize1,'Deflate',9);
            h5writeatt(h5OutputFile,[groupName,fieldNames{iField}],...
                'Dimensions','nImageSize x nImageSize');
        end
        
        iField = find(strcmp(fieldNames,'time'));
        h5create(h5OutputFile,[groupName,fieldNames{iField}],...
            [Inf Inf],'ChunkSize',chunkSize5,'Deflate',9);
        h5writeatt(h5OutputFile,[groupName,fieldNames{iField}],...
            'Dimensions','nTime x 1');
        h5writeatt(h5OutputFile,[groupName,fieldNames{iField}],...
            'Units','Posix time');
        [~, datasetInfo]=ish5dataset(h5OutputFile,groupName);
    end
    
    datasetNames=struct2cell(datasetInfo.Datasets);
    
    % Write ASI Images
    iField=find(strcmp(datasetNames(1,:),'ASI'));
    start = start + datasetInfo.Datasets(iField).Dataspace.Size;
    h5write(h5OutputFile,[groupName,'ASI'],...
                    ASI,start,count);
    % Write time
    iField=find(strcmp(datasetNames(1,:),'time'));
    start5 = start5 + datasetInfo.Datasets(iField).Dataspace.Size;
    h5write(h5OutputFile,[groupName,datasetNames{1,iField}],...
                        (input.time)',start5,count5);
    
   % Write coordinates, if they don't already exists
   [~,fields2D] =  intersect(fieldNames,{'lat','lon','az','el',...
            'mlat','mlon','alt'});
    for iField = fields2D'
        if ~ish5dataset(h5OutputFile,[groupName,fieldNames{iField}])
        h5write(h5OutputFile,[groupName,fieldNames{iField}],...
                (input.(fieldNames{iField})));  
        end
    end 
    
    % Write sensor location, if they already don't exist
    if ~ish5dataset(h5OutputFile,[groupName,'sensorloc'])
    dset = input.sensorloc;
    dset_details.Name = 'sensorloc';
    dset_details.Location = groupName;
    attr = '3x1 [lat,lon,alt-meters]';
    attr_details.Name = 'Dimensions';
    attr_details.AttachedTo = [groupName,dset_details.Name];
    attr_details.AttachType = 'dataset';
    hdf5write(h5OutputFile,dset_details,dset,...
        attr_details,attr,'WriteMode','append');
    end
        
    % Write/over-write creation date (if they already don't exist)/ last-updated duration
    startTime = unixtime2matlab(input.time(1));
    endTime = unixtime2matlab(input.time(end));
    
    if ~datasetExists
    h5writeatt(h5OutputFile,groupName,...
        'creation_date',datestr(now));
    else
    h5writeatt(h5OutputFile,groupName,...
        'updated_date',datestr(now));
    end
    h5writeatt(h5OutputFile,groupName,...
        'last_updated_duration_in_file',[datestr(startTime),' - ',datestr(endTime)]);    
    
end
