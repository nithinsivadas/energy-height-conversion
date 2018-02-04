function [ varargout ]=readMadhdf5(filename, paramstr)
%% readMad_hdf5
% by John Swoboda and Anna Stuhlmacher
% madgrigal h5 read in function for the python implementation of GeoData for Madrigal Sondrestrom data
% Input
% filename: path to hdf5 file
% paramstr: list of parameters to look at written as strings
% Returns:
% dictionary with keys are the Madrigal parameter string, the value is an array
% rows are unique data locations (data_loc) = (rng, azm, el1)
% columns are unique times
% 


%open hdf5 file
all_data = h5read(filename,'/Data/Table Layout');
sensor_struct = h5read(filename,'/Metadata/Experiment Parameters');
sensor_data = sensor_struct.value';

sensorname = sensor_data(1,:);

if strmatch('Sondrestrom', sensorname)
    radar = 1;
    disp('Sondrestrom data')
elseif strmatch('Poker Flat' , sensorname)
    radar = 2;
    disp('PFISR data')
elseif strmatch('Resolute Bay' , sensorname)
    radar = 3;
    disp('RISR data')
else
    error('Sensor type not supported by program in this version')   
end
%get the data location (range, el1, azm)

if radar == 1
    angle1 = 'elm';
    rng = all_data.('gdalt');
elseif radar == 2
    angle1 = 'elm';
    rng = all_data.('range');
elseif radar ==3
    angle1 = 'elm';
    rng = all_data.('range');
end

try
    el = all_data.(angle1);
catch 
    el = NaN(size(rng));
end

try
    azm = all_data.('azm');
catch 
    azm = NaN(size(rng));
end
% take out nans
nan_ar = isnan(rng)|isnan(el)|isnan(azm);
notnan = ~nan_ar;
all_loc= zeros(sum(notnan),3);

icount=1;
for i =1:length(rng)
    if notnan(i)
        all_loc(icount,:) = [rng(i),azm(i),el(i)];
        icount=icount+1;
    end
end

%create list of unique data location lists
[dataloc,~,icloc] = unique(all_loc,'rows');
times1 = all_data.('ut1_unix')(notnan);
times2 = all_data.('ut2_unix')(notnan);
all_times = [times1,times2];

[uniq_times,~,ictime] = unique(all_times,'rows');

%initialize and fill data dictionary with parameter arrays
data = struct();
maxcols = size(uniq_times,1);
maxrows = size(dataloc,1);
for ip =1:length( paramstr)
    p=paramstr{ip};
    if isempty(strmatch(p, fieldnames(all_data)))
        warning( [ p,  ' is not a valid parameter name.'])
        continue
    end
    tempdata = all_data.(p)(notnan); %list of parameter pulled from all_data
    temparray = zeros([maxrows,maxcols]); %converting the tempdata list into array form
    linind = sub2ind(size(temparray),icloc,ictime);
    temparray(linind)=double(tempdata);
    
    data.(p)=temparray;
end
%get the sensor location (lat, long, rng)
lat = str2double(sensor_data(8,:));
lon = str2double(sensor_data(9,:));
sensor_alt = str2double(sensor_data(10,:));
sensorloc =[lat,lon,sensor_alt];
coordnames = 'Spherical';

varargout = {data,coordnames,dataloc,sensorloc,double(uniq_times)};