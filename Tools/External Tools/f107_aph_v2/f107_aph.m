function [F107A, F107, APH] = f107_aph(daten,setForceDownload)

% INPUTS
% daten      : m x 1 array of datenum
% setForceDownload : Boolean true or false option to force download of data
%                    from amisr website

% OUTPUTS
% 'F107A'     : m x 1 array of 81-day average solar flux
% 'F107'      : m x 1 array of daily solar flux
% 'APH'       : m x 7 array of magnetic indices

% Computes/retrieves the 7 magnetic indices and two solar flux values
% (previous day and centered 81-day mean) for atmosnrlmsise00

% Author : John A. Smith, CIRES/NOAA, Univ. of Colorado at Boulder

% convert year, doy and UTseconds to column vectors
% year=reshape(year,[],1);
% doy=reshape(doy,[],1);
% UTseconds=reshape(UTseconds,[],1);

if nargin<2
    setForceDownload=false;
end

daten = reshape(daten,[],1);
UTseconds=seconds(timeofday(datetime(datevec(datenum(daten)))));

 
M = length(daten);

% acknowledge file ID #24235 "doy2date.m" for this insert
z = zeros(M,5);

% check that date is not in the future
if any(daten>now)
    error('Invalid time or date');
end

% determine year for 3 days prior (Ap)
[YAp, ~, ~, ~, ~, ~] = datevec(daten-3);

% determine year for 40 days prior
[YSF, ~, ~, ~, ~, ~] = datevec(daten-40);

% check if both Ap and solar flux (SF) data available for year
if min([YAp YSF])<1947
    error('Solar flux data only available from 1947 onward');
end

% convert year to string
yrstr = num2str(unique([YAp year(daten)]));
N = size(yrstr,1);

% define remote directories
url = 'https://amisr.com/geophys_params/';
remote = cellstr(strcat(url,yrstr));

% define full local directories
localFolder = strcat(initialize_root_path,'LargeFiles',filesep,'KP_AP',filesep);
local = cellstr(strcat(localFolder,yrstr));

% determine if needed Ap data already exists in pwd, if not then download
options = weboptions;
options.CertificateFilename=('');
for i=1:N
    if ~exist(local{i},'file') || setForceDownload
        websave(local{i},remote{i},options);% download dataset
    else
        warning(['File ',local{i},' alread exists']);
    end
end


X=[];
Y=[];
Z=[];

for j=1:N
    [x, y, z] = read_magnetic(local{j});
    X=[X;x];
    Y=[Y;y];
    Z=[Z;z];
end

APH=-1*ones(M,7); % preallocate APH matrix

% define and construct APH matrix
if (X(end)<daten)
    error(['The file ',url,datestr(daten,'yyyy'),' has only entries up to ',datestr(X(end))]);
end

for k=1:M
    row=find(X==floor(daten(k)));
    sub=Y(row-3:row,:)';
    ti=ceil(UTseconds(k)/10800+.001);
    APH(k,:)=[Z(row) sub(24+ti) sub(23+ti) sub(22+ti) ... 
        sub(21+ti) nanmean(sub((13:20)+ti)) nanmean(sub((5:12)+ti))];
end

% ----- finished magnetic index, starting solar flux -----

% create tracking file if not present (first time running)
if ~exist([localFolder,'f107.txt'],'file')
    fid = fopen([localFolder,'f107.txt'],'a');
    fwrite(fid,'0')
    fclose(fid);
end

% read date of previous file's modification
fid = fopen([localFolder,'f107.txt']);update = fscanf(fid,'%f');fclose(fid);

% determine when FTP server file last modified
f=ftp('ftp.ngdc.noaa.gov');
solar_dir = '/STP/space-weather/solar-data/solar-features/solar-radio/noontime-flux/penticton/penticton_observed/listings/';
cd(f,solar_dir);
latest = dir(f);
close(f);
last = latest(2).datenum;
filename = [localFolder 'listing_drao_noontime-flux-observed_daily.txt'];

% update if new data available
if ~exist(filename,'file') || last>(update+1000)
    f=ftp('ftp.ngdc.noaa.gov'); % open ftp session
    disp('downloading solar data...') 
    cd(f,solar_dir);
    mget(f,'listing_drao_noontime-flux-observed_daily.txt',localFolder);
    close(f); % close session
    fid = fopen([localFolder,'f107.txt'],'w');fprintf(fid,'%f',last);fclose(fid);
end

[x, y] = read_solarflux(filename); % read solar flux data

if x(end)<(max(daten)+40)
    warning('on')
    warning('81-day centered average not possible using available data')
end

warning('off') % turn off NaN warnings during interp

F107 = interp1(x,y,daten-1,'spline'); % interp F10.7 for requested dates

% calculate F10.7 81-day centered mean about doy
F107A = -1*ones(M,1);
for i=1:M
    index = x>(daten(i)-40) & x<(daten(i)+40);
    F107A(i) = nanmean(y(index));
end

end