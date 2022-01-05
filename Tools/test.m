    dkTime = 10;
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
            ASI(k,:,:) = fitsread(pathStr(iTime));
            time(k,1) = posixtime(datetime(fitsfiletimestamp(localFileListName(iTime)),'ConvertFrom','datenum'));
            k = k+1;
        end
        write_h5_dataset(h5FileStr,[datasetPath,'time'],time,1,true);
        write_h5_dataset(h5FileStr,[datasetPath,'ASI'],ASI,1,true);
    end