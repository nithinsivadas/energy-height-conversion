function read_all_local_FITS_and_create_HDF5(h5FileStr, wavelengthStr, tempStorePath)
    localFileList = dir(strcat(tempStorePath,filesep,'*.FITS'));
    localFileListName = string(deblank(char(localFileList.name)));
    pathStr = string(strcat(deblank(char(localFileList.folder)),filesep,localFileListName));
    
    dkTime = 800;
    nTimeASITotal = length(pathStr);
    nkTime = ceil(nTimeASITotal/dkTime);
    
    datasetPath = char(strcat('/',wavelengthStr,'/'));
    minChunkTimeDim = min(mod(nTimeASITotal,dkTime),10);
    for kTime=1:1:nkTime
        timeEndIndx = min(kTime*dkTime,nTimeASITotal);
        timeStartIndx = 1 + (kTime-1)*dkTime;
        k=1;
        ASI = [];
        time = [];
        errorFlag = [];
        for iTime = timeStartIndx:1:timeEndIndx
            try
                ASI(k,:,:) = fitsread(pathStr(iTime));
                errorFlag(k,1) = 0;
            catch ME
                ASI(k,:,:) = nan;
                errorFlag(k,1) = 1;
            end
            try
                time(k,1) = posixtime(datetime(fitsfiletimestamp(localFileListName(iTime)),'ConvertFrom','datenum'));
                errorFlag(k,1) = 0;
            catch ME
                time(k,1) = nan;
                errorFlag(k,1) = 2;
            end
            k = k+1;
        end      
        write_h5_dataset(h5FileStr,[datasetPath,'time'],time,1,true,minChunkTimeDim,false,9);
        write_h5_dataset(h5FileStr,[datasetPath,'ASI'],ASI,1,true,minChunkTimeDim,false,5);
        write_h5_dataset(h5FileStr,[datasetPath,'errorFlag'],errorFlag,1,true,minChunkTimeDim,false,9);
    end
    
end