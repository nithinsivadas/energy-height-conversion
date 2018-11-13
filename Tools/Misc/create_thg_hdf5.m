function status = create_thg_hdf5(siteName,outputH5FileStr,varargin)
%create_thg_hdf5.m Downloads, and writes themis GBO ASI images for a
%particular site to and HDF5 file. If time limits are set as default, then
%the function will search the outputH5FileStr for a time array and use that
%to download the right cdf files from themis webserver. 

p = inputParser;
addParameter(p,'dateFormat','yyyy-mm-dd HH:MM:SS',@(x) isstring(x)||ischar(x));
addParameter(p,'minTimeStr','default',@(x) isstring(x)||ischar(x));
addParameter(p,'maxTimeStr','default',@(x) isstring(x)||ischar(x));
addParameter(p,'localStorePath','default',@(x) isstring(x)||ischar(x));
addParameter(p,'sensorLocationFile','default',@(x) isstring(x)||ischar(x));
addParameter(p,'altIndx',2,@(x) x>=1&&x<=3&&(ceil(x)-x)==0);

addRequired(p,'outputH5FileStr',@(x)contains(x,{'.h5','.hdf5'})); 
addRequired(p,'siteName', @(x) isstring(x)||ischar(x));

parse(p,outputH5FileStr,siteName,varargin{:});

%% Initialization
minTimeStr = p.Results.minTimeStr;
if strcmp(minTimeStr,'default')
    if ish5dataset(outputH5FileStr,'/energyFluxFromMaxEnt/time')
        time = h5read(outputH5FileStr,'/energyFluxFromMaxEnt/time');
        minTimeStr = datestr(time(1),p.Results.dateFormat);
    else
        error('23: h5 file does not have time /energyFluxFromMaxEnt/time dataset \n must input minTimeStr');
    end
end

maxTimeStr = p.Results.maxTimeStr;
if strcmp(maxTimeStr,'default')
    if ~exist('time')
        if ish5dataset(outputH5FileStr,'/energyFluxFromMaxEnt/time')
            time = h5read(outputH5FileStr,'/energyFluxFromMaxEnt/time');
            maxTimeStr = datestr(time(1),p.Results.dateFormat);
        else
            error('33: h5 file does not have time /energyFluxFromMaxEnt/time dataset \n must input maxTimeStr');
        end
    else
        maxTimeStr = datestr(time(end),p.Results.dateFormat);
    end
end

sensorLocationFile = p.Results.sensorLocationFile;
if strcmp(sensorLocationFile,'default')
    sensorLocationFile = [initialize_root_path,'energy-height-conversion',...
        filesep,'Tools',filesep,'External Tools',filesep,'thmasi',filesep,...
        'THEMIS_ASI_Station_List_Nov_2011.xls'];
end

%% Download the data
[cdfFileList, cdfCalFile, status, cmdout] = ...
    download_thg(minTimeStr, maxTimeStr, siteName,...
    'dateFormat',p.Results.dateFormat);

nFiles = length(cdfFileList);
thgSites =  parse_thg_location_xls(sensorLocationFile);
siteID = find(strcmpi(thgSites.code,siteName));
if isempty(siteID)
    error(['Site - ',siteName,' not in the sensorLocationFile: ',...
        sensorLocationFile]);
end
multiWaitbar('Writing Themis GBO ASI to HDF5 Files',0);
id = 1./nFiles;
for iFile = 1:1:nFiles
    multiWaitbar('Writing Themis GBO ASI to HDF5 Files','Increment',id);
    thgData = parse_thg_cdfData(char(cdfFileList(iFile)),char(cdfCalFile));
    if(iFile == 1)
        disp(['Choosing altitude to be ',...
            num2str(thgData.alt(p.Results.altIndx)/1000),' km']);
    end
    write_thg_to_hdf5(outputH5FileStr,thgData.ASI,'lat',thgData.glat,...
    'lon',thgData.glon,'az',thgData.az,'el',thgData.el,'alt',thgData.alt,...
    'altIndx',p.Results.altIndx,'mlat',thgData.mlat,'mlon',thgData.mlon,...
    'time',thgData.time,'sensorloc',[thgSites.glat(siteID),thgSites.glon(siteID),0],...
    'siteCode',siteName);        
end
add_ASI_background_to_hdf5(siteName,outputH5FileStr); % Adding the backgroud
multiWaitbar('Writing Themis GBO ASI to HDF5 Files',1);
status = 'complete';    
end


