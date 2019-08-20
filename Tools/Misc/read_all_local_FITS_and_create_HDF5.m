function read_all_local_FITS_and_create_HDF5(h5FileStr, wavelengthStr, folderStr, tempStorePath)
    localFileList = dir(strcat(tempStorePath,filesep,'*.FITS'));
    localFileListName = string(deblank(char(localFileList.name)));
    pathStr = string(strcat(deblank(char(localFileList.folder)),filesep,localFileListName));
    
    dkTime = 800;
    nTimeASITotal = length(pathStr);
    nkTime = ceil(nTimeASITotal/dkTime);
    
    datasetPath = char(strcat('/',folderStr,'/',wavelengthStr,'/'));
    for kTime=1:1:nkTime
        timeEndIndx = min(kTime*dkTime,nTimeASITotal);
        timeStartIndx = 1 + (kTime-1)*dkTime;
        k=1;
        ASI = [];
        time = [];
        for iTime = timeStartIndx:1:timeEndIndx
            try
            ASI(k,:,:) = fitsread(pathStr(iTime));
            catch ME
            ASI(k,:,:) = nan;
            end
            time(k,1) = posixtime(datetime(fitsfiletimestamp(localFileListName(iTime)),'ConvertFrom','datenum'));
            k = k+1;
        end      
        write_h5_dataset(h5FileStr,[datasetPath,'time'],time,1,true);
        write_h5_dataset(h5FileStr,[datasetPath,'ASI'],ASI,1,true);
    end
    
end