function [cdfFileList, cdfCalFile, status, cmdout] = download_thg(minTimeStr, maxTimeStr, siteName, varargin)
    % download_thg.m Downloads themis GBO all sky imager data (asf) to
    % local folder, without overlap with files already in folder. 
    
    p = inputParser;
    
    addParameter(p,'localStorePath','default',@(x) isstring(x)||ischar(x));
    addParameter(p,'dateFormat','yyyy-mm-dd HH:MM:SS',@(x) isstring(x)||ischar(x));
    
    addRequired(p, 'minTimeStr', @(x) isstring(x)||ischar(x));
    addRequired(p, 'maxTimeStr', @(x) isstring(x)||ischar(x));
    addRequired(p,'siteName', @(x) isstring(x)||ischar(x));
    parse(p,minTimeStr,maxTimeStr,siteName,varargin{:});
    
    if strcmp(p.Results.localStorePath,'default')
        localStorePath = [initialize_root_path,'LargeFiles',filesep,'thg',filesep];
    else
        localStorePath = p.Results.localStoreFolder;
    end
    
    
    host = 'http://themis.ssl.berkeley.edu';
    remoteLink = '/data/themis/thg/l1/asi/';
    arrayYYMM = parse_yy_mm(minTimeStr, maxTimeStr, 'Format', p.Results.dateFormat);
    fileNames = get_thg_fileNames(minTimeStr, maxTimeStr, p.Results.dateFormat, siteName)';
    thgFileNames = fileNames.thgFileNames;
    thgCalFile = fileNames.thgCalFile;
    
    % Check if local store Path exists
    if(isdir(localStorePath))
        localFileList = dir(localStorePath);
        localFileListName = strtrim(string(char(localFileList.name)));
    else
        mkdir(localStorePath);
        localFileListName = ".";
    end
    
    % Generate all the download links
    remoteFileList = [];
    remoteFileListLink = [];
    for i=1:length(arrayYYMM(2,:))
        remoteLinkFinal = [remoteLink,siteName,'/',num2str(arrayYYMM(1,i)),'/',num2str(arrayYYMM(2,i),'%02d')];
        S = webread([host,remoteLinkFinal],'option','text');
        remoteFileListTemp = extractBetween(S,' ]"></td><td><a href="','">thg');
        if ~isempty(remoteFileListTemp)
        remoteFileList = [remoteFileList; remoteFileListTemp];
        remoteFileListLink = [remoteFileListLink; strcat(string(remoteLinkFinal),'/',string(char(remoteFileListTemp)))]; 
        end
    end
    
    % Generate list of files already downloaded
    localFileList = dir(localStorePath);
    localFileListName = strtrim(string(char(localFileList.name)));
    
    try %(If there are no files online then there will be an error)
        % Identifying the asf files to be downloaded
        [intersectFiles,toDownloadIndx] = intersect(remoteFileList,thgFileNames);
        intersectFileLinks = remoteFileListLink(toDownloadIndx);
        cdfFileList = strcat(localStorePath,intersectFiles);
        
        % Add Calibration File
        cdfCalFile = strcat(localStorePath, thgCalFile);
        intersectFileLinks = [intersectFileLinks; strcat("/data/themis/thg/l2/asi/cal/",thgCalFile)];
        intersectFiles = [intersectFiles; thgCalFile];

        % Identifying identical files in remote and local directories
        [~, ia] = setdiff(intersectFiles,localFileListName);
        toDownloadFileLinks = intersectFileLinks(ia);

        urlFile = 'tempURL.txt';
        urlFilePath = [localStorePath,filesep,urlFile];
        urls = strcat(host,toDownloadFileLinks);
        fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
        if isunix
        [status,cmdout]=unix(['aria2c -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
        else
        [status,cmdout]=system(['aria2c -V -c -j 50 ','-d ',localStorePath,' -i ',urlFilePath]);
        end
    
    catch ME
        error('Possibly no files online during this time period.');
    end
end

function [data]=get_thg_fileNames (minTimeStr, maxTimeStr, format, siteName)
    % Generates the ideal file name (l1, asf, v01)
    startyy = year(minTimeStr, format);
    startmm = month(minTimeStr, format);
    startdd = day(minTimeStr,format);
    starthh = hour(minTimeStr,format);
    
    endyy = year(maxTimeStr, format);
    endmm = month(maxTimeStr, format);
    enddd = day(maxTimeStr,format);
    endhh = hour(maxTimeStr,format);
    
    k = 1;
    ihh=starthh;
    while datetime(startyy,startmm,startdd,ihh,0,0)<=datetime(endyy,endmm,enddd,endhh,0,0)
        datetimeStr=datestr(datetime(startyy,startmm,startdd,ihh,0,0),'yyyymmddHH');
        thgFileNames(k)=string(['thg_l1_asf_',siteName,'_',datetimeStr,'_v01.cdf']);
        ihh = ihh+1;
        k = k+1;
    end
    thgCalFile = string(['thg_l2_asc_',siteName,'_19700101_v01.cdf']);
    data.thgFileNames=thgFileNames';
    data.thgCalFile = thgCalFile;
end