function [ varargout ] = readSRI_h5( varargin )
%% Description: read_SRI
% By Andrew Lee
% This function will read the SRI format data into the Geodata format structure.
% Input
% filename - A string for the filename with the data.
% parameter - A cell array with the parameters names that are desired.
% timebounds - A vector of start and end times in posix.
% Output
% Varargout - A cell array with the GeoData values. See GeoData
% Documentation for details.
% {'data','coordnames','dataloc','sensorloc','times'};
%% Variable Arguments Out
varnames = {'data','coordnames','dataloc','sensorloc','times'};
varargout = cell(1,length(varnames));
%% Variable Arguements In
filename = varargin{1};
parameter = varargin{2};


%% Load SRI H5 file

datadict = {'Ne', 'dNe', 'Vi', 'dVi', 'Ti', 'dTi', 'Te', 'dTe';
    '/FittedParams/Ne', '/FittedParams/dNe', '/FittedParams/Fits', '/FittedParams/Errors', '/FittedParams/Fits',...
    '/FittedParams/Errors', '/FittedParams/Fits', '/FittedParams/Errors'};

% Define coordinate type
coordnames = 'Spherical';
varargout{2} = coordnames;

%% Load times
loadtimes = double(hdf5read(filename,'/Time/UnixTime')');

if nargin <3
    times = loadtimes;
    T1 = 1;
    T2 = size(loadtimes,1);
else 
    timebound = varargin{3};
    if isa(timebound,'cell')
       timebound = datestr2unix(timebound);
    end
    T1 = find(loadtimes(:,1) >= timebound(1), 1, 'first');
    T2 = find(loadtimes(:,2) >= timebound(2), 1, 'first');
    times = loadtimes(T1:T2,:);
end

% Select and load desired times
varargout{5} = times;

%% Load sensor locations
lat = hdf5read(filename,'/Site/Latitude');
lon = hdf5read(filename,'/Site/Longitude');
alt = hdf5read(filename,'/Site/Altitude');
sensorloc = [lat lon alt];
varargout{4} = sensorloc;
%% Load Data Locations
range = hdf5read(filename,'/FittedParams/Range')/1000; %in km
angles = hdf5read(filename,'BeamCodes');
angles = angles(2:3,:);
nrange = length(range(:,1));
repangles = repmat(angles,1,nrange);
allaz = repangles(1:2:end);
allel = repangles(2:2:end);
tempdataloc = [reshape(range.',1,[]);allaz;allel];
dataloc = tempdataloc.';
varargout{3} = dataloc;

%% Load Data
% paramindex = find(strcmp(datadict(1,:),parameter)==1);
data = struct();
for k = 1:length(parameter)
    paramindex = find(strcmp(datadict(1,:),parameter{k}));
    
    if paramindex == 0
         error('provide correct parameter (e.g. ''Ne'')')
    else
        tempData = hdf5read(filename,char(datadict(2,paramindex)));
    end
    
    if (paramindex == 3 || paramindex == 4) %For Vi or dVi
        tempData = squeeze(tempData(4,1,:,:,:));
    elseif (paramindex == 5 || paramindex == 6) %For Ti or dTi
        tempData = squeeze(tempData(2,1,:,:,:));
    elseif (paramindex == 7 || paramindex == 8) %For Te or dTe
        tempData = squeeze(tempData(2,end,:,:,:));
    end
    tempData = tempData(:,:,T1:T2);
    data.(parameter{k}) = reshape(permute(tempData,[2 1 3]),[size(dataloc,1),size(times,1)]);
end
varargout{1} = data;