function data = add_temp_noaa17_hdf5(inputMatFile, outputH5FileStr,omniH5FileStr, varargin)
%add_temp_noaa17_hdf5 Adds particular noaa17 mat file, also creates required
% magnetic foot points specified by magneticFieldModel parameter. 
% This is a temporary, function, and can be easily modified to download,
% run fortran code, and automatically get POES data. 

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};

addParameter(p,'dateFormat','yyyy-mm-dd HH:MM:SS',@(x) isstring(x)||ischar(x));
addParameter(p,'minTimeStr','default',@(x) isstring(x)||ischar(x));
addParameter(p,'maxTimeStr','default',@(x) isstring(x)||ischar(x));
addParameter(p,'magneticFieldModel','TS89',@(x) any(strcmpi(x,expectedMagFieldModels)));
addParameter(p,'localStorePath','default',@(x) isstring(x)||ischar(x));
addParameter(p,'projectionAltitude',110);

addRequired(p,'inputMatFile', @(x)contains(x,{'.mat'}));
addRequired(p,'outputH5FileStr',@(x)contains(x,{'.h5','.hdf5'})); 
addRequired(p,'omniH5FileStr',@(x)contains(x,{'.h5','.hdf5'})); 

parse(p,inputMatFile,outputH5FileStr,omniH5FileStr,varargin{:});

%% Initialization
minTimeStr = p.Results.minTimeStr;
if strcmp(minTimeStr,'default')
    if ish5dataset(outputH5FileStr,'/energyFluxFromMaxEnt/time')
        time = h5read(outputH5FileStr,'/energyFluxFromMaxEnt/time');
        minTimeStr = datestr(time(1),p.Results.dateFormat);
    else
        error('23: h5 file does not have time /energyFluxFromMaxEnt/time dataset. Must input minTimeStr');
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

%% Download NOAA data
load(inputMatFile);
[maginput,timeMaginput]=generate_maginput(omniH5FileStr,minTimeStr,maxTimeStr);

poesTime = poes.time;
poesTimeIndx = find_time(poesTime,minTimeStr):1:find_time(poesTime,maxTimeStr);

data.time = poesTime(poesTimeIndx);
data.GDZ = [poes.lat(poesTimeIndx,1), poes.lon(poesTimeIndx,1), 850*ones(size(poes.lat(poesTimeIndx,1)))]; 
data.XYZ_GEO = onera_desp_lib_rotate([data.GDZ(:,3)...
    ,data.GDZ(:,1), data.GDZ(:,2)],'gdz2geo');
data.magFieldStr = p.Results.magneticFieldModel;
magFieldNo = find_irbem_magFieldModelNo(data.magFieldStr);
sysaxes = 1; % GEO

nTime = length(data.time);
for iTime = 1:1:nTime
    thisMaginput = interp1(timeMaginput',maginput,data.time(iTime));
    data.NFoot(iTime,:) = geopack_find_foot_point...
            (magFieldNo,100,sysaxes,data.time(iTime),...
            data.XYZ_GEO(iTime,1),...
            data.XYZ_GEO(iTime,2),...
            data.XYZ_GEO(iTime,3),...
            p.Results.projectionAltitude,+1,thisMaginput); % IRBEM Has issues converging?
     data.Lm(iTime,:) = onera_desp_lib_make_lstar...
            (magFieldNo,[0,0,0,0,0],sysaxes,data.time(iTime),...
            data.XYZ_GEO(iTime,1),...
            data.XYZ_GEO(iTime,2),...
            data.XYZ_GEO(iTime,3),...
            thisMaginput);   
end

write_sc_to_hdf5(outputH5FileStr,'noaa17',data.time,data.XYZ_GEO,...
    data.magFieldStr,data.NFoot,data.Lm);

end

