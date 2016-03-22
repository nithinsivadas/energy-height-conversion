function [density, temperature, f_10_7_used, Ap_used] = msis(time, ...
    latitude, longitude, altitude, utc, coord, curldir, f_10_7_daily, ...
    f_10_7_3month, ap_daily)

% MSIS Atmosphere model MSIS-E-90.
% 
% Usage: [DENSITY, TEMPERATURE, F10_7_USED, AP_USED] = MSIS(TIME, LATITUDE,
%             LONGITUDE, ALTITUDE, UTC, COORD, CURLDIR, F10_7_DAILY,
%             F10_7_3MONTH, AP_DAILY)
% 
% Computes the MSIS-E-90 Atmosphere Model, which is an model for Earth's
% atmosphere from the ground into the exosphere. MSIS stands for Mass
% Spectrometer and Incoherent Scatter Radar, which are the two data sources
% used to develop the model. The position and time inputs can be scalars or
% arrays; when they are arrays, they should all have the same number of
% elements. The outputs will be arrays with the same number of rows as
% elements in the input arrays, possibly 1 (so the shape of the input is
% not preserved).
% 
% The function makes the computation by querying the online interface at
% http://omniweb.gsfc.nasa.gov/vitmo/msis_vitmo.html (hence internet access
% is required), which makes it pretty slow, especially when either more
% than one of the position and time inputs are made to vary or if the
% position and time inputs are not spaced linearly. If more than one input
% is varied, the function can be sped up by using a for loop and holding
% the smaller arrays constant, assuming the largest array is spaced
% linearly. This will result in fewer calls to the website since the
% website allows for a linear sweep in one variable. See the script
% msistest.m for an example of this. There is one exception: The online
% interface has odd behavior that does not allow for sweeps in longitude
% for any altitude except the default (100 km), so longitude sweeps will be
% computed one at a time unless ALTITUDE is 100.
% 
% The query is made using the command curl in an operating system terminal.
% This is built-in to Unix but not Windows. curl for Windows can be
% downloaded from http://curl.haxx.se/download.html. The directory where
% the curl.exe file can be found should be passed into CURLDIR for Windows
% computers. CURLDIR defaults to the same directory as this function.
% 
% A value for -1 means the output is invalid for the given input.
% 
% This is NOT the most recent MSIS model. The website for that version is:
% http://www.nrl.navy.mil/research/nrl-review/2003/atmospheric-science/picone/
% 
% Inputs:
%   -TIME: Times to compute MSIS model either in MATLAB serial date number
%   format or a string that can be converted into MATLAB serial date number
%   format using DATENUM with no format specified (see documentation of
%   DATENUM for more information). Whether the times are local or UTC are
%   determined by the input UTC. Valid range is from year 1958 to year 2013
%   currently (optional, default is January 1, 2000 at 01:30).
%   -LATITUDE: Latitude in degrees to compute MSIS model. Whether this is
%   geodetic, geocentric, or geomagnetic latitude is determined by the
%   input COORD. Valid range is -90 degrees to 90 degrees (optional,
%   default is 50 degrees).
%   -LONGITUDE: Longitude in degrees to compute MSIS model. Whether this is
%   geodetic, geocentric, or geomagnetic longitude is determined by the
%   input COORD. Valid range is 0 degrees to 360 degrees (optional, default
%   is 40 degrees).
%   -ALTITUDE: For geodetic or geomagnetic coordiates, the height in km
%   above the Earth's surface. For geocentric coordiates, the radius in km
%   from the center of the Earth. Valid range for altitude is 0 km to 1000
%   km (optional, default when all other inputs are scalars is to sweep
%   from 0:50:1000 km and when any other input is an array is 100 km).
%   -UTC: Set to true (or 'UTC', 'U') if the times in TIME are in
%   Coordinated Universal Time (UTC) or false (or 'Local', 'LT', 'L') if
%   the times in TIME are in local time (optional, default is true).
%   -COORD: String specifying the coordinate system to use. Can be
%   geodetic ('geodetic', 'geod', or 'gd'), geomagnetic ('geomagnetic',
%   'geom', or 'gm'), or geocentric ('geocentric', 'geom', or 'gm')
%   (optional, default is geodetic).
%   -CURLDIR: Directory where curl.exe can be found (optional, only
%   necessary for Windows computers, and default for those is the same
%   directory that this function is located).
%   -F_10_7_DAILY: F_10.7 daily index to use in the model. Valid range is 0
%   to 400 (optional, default is to leave this field blank, in which case
%   it is taken from "real data base" [sic]).
%   -F_10_7_3MONTH: F_10.7 3 month average index to use in the model. Valid
%   range is 50 to 350 (optional, default is to leave this field blank, in
%   which case it is taken from "real data base" [sic]).
%   -AP_DAILY: Daily Ap index to use in the model. Valid range is 0 to 40
%   (optional, default is to leave this field blank, in which case it is
%   taken from "real data base" [sic]).
% 
% Outputs:
%   -DENSITY: Array with each column having the following:
%       1. O: Atomic oxygen (O) number density in m^-3.
%       2. N2: Molecular nitrogen (N_2) number density in m^-3.
%       3. O2: Molecular oxygen (O_2) number density in m^-3.
%       4. MASS_DENSITY: Total mas density in kg/m^3.
%       5. HE: Atomic helium (He) number density in m^-3.
%       6. AR: Atomic argon (Ar) number density in m^-3.
%       7. H: Atomic hydrogen (H) number density in m^-3.
%       8. N: Atomic nitrogen (N) number density in m^-3.
%   -TEMPERATURE: Array with each column having the following:
%       1. TN: Neutrals temperature in K.
%       2. TEX: Exospheric temperature in K.
%   -F10_7_USED: Array with each column having the following:
%       1. F10.7 daily index used in the model.
%       2. F10.7 3 month average index used in the model.
%   -AP_USED: Array with each column having the following:
%       1. Daily Ap index used in the model.
%       2. Ap index used in the model from 0 to 3 hours prior.
%       3. Ap index used in the model from 3 to 6 hours prior.
%       4. Ap index used in the model from 6 to 9 hours prior.
%       5. Ap index used in the model from 9 to 12 hours prior.
%       6. Ap index used in the model from 12 to 33 hours prior.
%       7. Ap index used in the model from 33 to 59 hours prior.
% 
% See also: MSISTEST, IRI, IGRF, DATENUM.

% Directory of this function.
fpath = mfilename('fullpath');
fpath = fpath(1:end-length(mfilename));

% Default behavior.
if nargin < 1 || isempty(time)
    time = datenum([2000 1 1 1 30 0]);
end
if nargin < 2 || isempty(latitude)
    latitude = 55;
end
if nargin < 3 || isempty(longitude)
    longitude = 45;
end
if nargin < 4 || isempty(altitude)
    if (ischar(time) || numel(time) == 1) && ...
            numel(latitude) == 1 && numel(longitude) == 1
        altitude = 0:50:1000;
    else
        altitude = 100;
    end
end
if nargin < 5 || isempty(utc)
    utc = true;
elseif ischar(utc)
    switch lower(utc)
        case {'utc', 'u'}
            utc = true;
        case {'local', 'lt', 'l'}
            utc = false;
        otherwise
            error('msis:badUTC', ['Unrecognized command ' utc '. Valid' ...
                ' options are ''utc'', ''u'', ''local'', ''lt'', or ' ...
                '''l''.']);
    end
end
if nargin < 6 || isempty(coord)
    coord = 'geod';
end
if isunix || ismac  % Curl is built-in to operating system.
    curldir = [];
elseif nargin < 7 || isempty(curldir)
    curldir = fpath;
end
if nargin < 8
    f_10_7_daily = [];
end
if nargin < 9
    f_10_7_3month = [];
end
if nargin < 10
    ap_daily = [];
end

% Convert time to a datenumber if it is a string.
if ischar(time)
    time = datenum(time);
end

% Convert the coordinates.
switch lower(coord)
    case {'geodetic', 'geod', 'gd'}
        geo_flag = '0.';
    case {'geocentric', 'geoc', 'gc'}
        % Convert coordinates to geodetic. The function ecef2geod assumes
        % meters, but we want km here.
        [x, y, z] = sph2cart(longitude*pi/180, latitude*pi/180, ...
            altitude*1e3);
        [latitude, longitude, altitude] = ecef2geod(x, y, z);
        altitude = altitude/1e3;
        geo_flag = '0.';
    case {'geomagnetic', 'geom', 'gm'}
        geo_flag = '1.';
    otherwise
        error('iri:coordCommandUnknown', ['Command ' coord ' unknown. ' ...
            'Valid options are ''geodetic'', ''geocentric'', and ' ...
            '''geomagnetic''.']);
end

% Error checking and input coversion.
if any(latitude < -90) || any(latitude > 90)
    error('msis:invalidLatitude', ['Input LATITUDE must be between ' ...
        '-90 degrees and 90 degrees.']);
end
longitude = mod(longitude, 360);
if any(altitude < 0) || any(altitude > 1e3)
    error('msis:invalidAltitude', ['Input ALTITUDE must be between ' ...
        '0 km and 1000 km.']);
end
if isempty(f_10_7_daily) || (f_10_7_daily >= 0 && f_10_7_daily <= 400)
    f_10_7_daily = sprintf('&f10_7=%#g', f_10_7_daily);
else
    error('msis:invalidF10_7_daily', ['Input F10_7_DAILY is %g but ' ...
        'must be between 0 and 400.'], f_10_7_daily);
end
if isempty(f_10_7_3month) || (f_10_7_3month >= 50 && f_10_7_3month <= 350)
    f_10_7_3month = sprintf('&f10_7_3=%#g', f_10_7_3month);
else
    error('msis:invalidF10_7_3month', ['Input F10_7_3MONTH is %g but ' ...
        'must be between 50 and 350.'], f_10_7_3month);
end
if isempty(ap_daily) || (ap_daily >= 0 && ap_daily <= 40)
    ap_daily = sprintf('&ap=%#g', ap_daily);
else
    error('iri:invalidTec_hmax', ['Input TEC_HMAX is %g km but must ' ...
        'be between 50 km and 2000 km.'], ap_daily);
end

% End of the string to be input into the system function.
endcmd = [f_10_7_daily f_10_7_3month ap_daily ...
    sprintf('&vars=%i', [5, 8:26]) ...
    '" http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi > "' ...
    ...'" http://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi > "' ...
    fullfile(fpath, 'temp.html"')];

% Get year, month, day, hour for MSIS.
[year, month, day, hour, minute, second] = datevec(time);
hour = hour + minute/60 + second/3600;
dayyear = dayofyear(year, month, day);

% Get the beginning of the string to be input into the system function.
if isunix || ismac % Curl is built-in to operating system.
    initialcmd = 'curl -d "model=msis';
else % Use the input curldir (possibly the default for that input).
    initialcmd = ['"' fullfile(curldir, 'curl"') ' -d "model=msis'];
end

% If one of the inputs can be swept keeping the others constant, run a
% sweep on the model.
A = [numel(unique(altitude)), numel(unique(latitude)), ...
    numel(unique(longitude)), numel(unique(year)), ...
    numel(unique(dayyear)), numel(unique(hour))];
if all(A == [1 1 1 1 1 1]) % Just one unique input.
    profile = 0;
elseif all(A == [A(1) 1 1 1 1 1]) % Sweep altitude.
    [sweep, tmp, sortind] = unique(altitude);
    if all(diff(diff(sweep)) < 1e-12)
        profile = 1;
        altitude = sweep(1);
        sweepmax = 1000;
    else % Altitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 A(2) 1 1 1 1]) % Sweep latitude.
    [sweep, tmp, sortind] = unique(latitude);
    if all(diff(diff(sweep)) < 1e-12)
        profile = 2;
        latitude = sweep(1);
        sweepmax = 90;
    else % Latitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 A(3) 1 1 1]) % Sweep longitude.
    [sweep, tmp, sortind] = unique(longitude);
    % Online interface can only have 100 km for longitude sweeps!
    if all(diff(diff(sweep)) < 1e-12) && altitude == 100
        profile = 3;
        longitude = sweep(3);
        sweepmax = 360;
    else % Longitude is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 A(4) 1 1]) % Sweep year.
    [sweep, tmp, sortind] = unique(year);
    if all(diff(diff(sweep)) < 1e-12)
        profile = 4;
        year = sweep(1);
        sweepmax = 2012;
    else % Year is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 1 A(5) 1]) % Sweep day of year or month.
    [sweep, tmp, sortind] = unique(dayyear);
    if all(diff(diff(sweep)) < 1e-12)
        profile = 7;
        dayyear = sweep(1);
        sweepmax = 366;
    % Day of year steps are not linear, but maybe month steps are.
    elseif numel(day(1)) == 1
        [sweep, tmp, sortind] = unique(month);
        profile = 5;
        month = sweep(1);
        sweepmax = 12;
    else % Day of year is unsweepable because steps are not linear.
        profile = 9;
    end
elseif all(A == [1 1 1 1 1 A(6)]) % Sweep hour in a day.
    [sweep, tmp, sortind] = unique(hour);
    if all(diff(diff(sweep)) < 1e-12)
        profile = 8;
        hour = sweep(1);
        sweepmax = 24;
    else % Hour is unsweepable because steps are not linear.
        profile = 9;
    end
else
    profile = 9;
end

% Call curl depending on the sweep profile.
switch profile
    case 0
        [status, result] = system([initialcmd ...
            sprintf('&year=%i', year) ...
            sprintf('&month=%i', month) ...
            sprintf('&day=%i', day) ...
            sprintf('&time_flag=%i', ~utc) ...
            sprintf('&hour=%#g', hour) ...
            '&geo_flag=' geo_flag ...
            sprintf('&latitude=%#g', latitude) ...
            sprintf('&longitude=%#g', longitude) ...
            sprintf('&height=%#g', altitude) ...
            '&profile=1' ...
            sprintf('&start=%#g', altitude) ...
            sprintf('&stop=%#g', altitude) ...
            '&step=1.' ...
            endcmd]);
        if status == 0
            data = parseresult(fpath, 1);
        else
            error('msis:curlError', ['Curl command did not work. ' ...
                'It returned status:%i, ' regexprep(result, '\\', ...
                '\\\\')], status);
        end
    case {1 2 3 4 5 7 8}
        % The online interface will only output up to a sweep length of
        % 500, so split the sweep into increments of 500.
        nsweeps = ceil(numel(sweep) / 500);
        sweepstart = sweep(1 : 500 : end);
        sweepstop = [sweep(500 : 500 : end), sweep(end)];
        sweepstep = mode(diff(sweep));
        sweeplen = round((sweepstop - sweepstart) ./ sweepstep + 1);
        prevsweeplen = [0 sweeplen(1:end-1)];
        sweepstop = min([sweep(end) + sweepstep/10, sweepmax]);
        data = zeros(20*sum(sweeplen), 1);
        for index = 1:nsweeps
            [status, result] = system([initialcmd ...
                sprintf('&year=%i', year(1)) ...
                sprintf('&month=%i', month(1)) ...
                sprintf('&day=%i', day(1)) ...
                sprintf('&time_flag=%i', ~utc) ...
                sprintf('&hour=%#g', hour(1)) ...
                '&geo_flag=' geo_flag ...
                sprintf('&latitude=%#g', latitude(1)) ...
                sprintf('&longitude=%#g', longitude(1)) ...
                sprintf('&height=%#g', altitude(1)) ...
                sprintf('&profile=%i', profile) ...
                sprintf('&start=%#g', sweepstart(index)) ...
                sprintf('&stop=%#g', sweepstop) ...
                sprintf('&step=%#g', sweepstep) ...
                endcmd]);
            if status == 0
                data((1:20*sweeplen(index)) + 20*prevsweeplen(index)*...
                    (index-1)) = parseresult(fpath, sweeplen(index));
            else
                error('msis:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, ' ...
                    regexprep(result, '\\', '\\\\')], status);
            end
        end
    % profile = 9 means there is more than one unique run to make. Turn all
    % the scalars into vectors with the same number of elements as the
    % largest array. If there is an array smaller than the largest array,
    % throw an error.
    case 9
        maxnum = max([numel(altitude), numel(latitude), ...
            numel(longitude), numel(year), numel(day), numel(month), ...
            numel(hour)]);
        if numel(altitude) == 1
            altitude = repmat(altitude, maxnum, 1);
        elseif numel(altitude) ~= maxnum && numel(altitude) > 1
            error('msis:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(latitude) == 1
            latitude = repmat(latitude, maxnum, 1);
        elseif numel(latitude) ~= maxnum && numel(latitude) > 1
            error('msis:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(longitude) == 1
            longitude = repmat(longitude, maxnum, 1);
        elseif numel(longitude) ~= maxnum && numel(longitude) > 1
            error('msis:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        if numel(time) == 1
            year = repmat(year, maxnum, 1);
            month = repmat(month, maxnum, 1);
            day = repmat(day, maxnum, 1);
            hour = repmat(hour, maxnum, 1);
        elseif numel(time) ~= maxnum && numel(time) > 1
            error('msis:invalidSize', ['Input vectors must all have ' ...
                'the same number of elements.']);
        end
        data = zeros(20*maxnum, 1);
        for index = 1:maxnum
            [status, result] = system([initialcmd ...
                sprintf('&year=%i', year(index)) ...
                sprintf('&month=%i', month(index)) ...
                sprintf('&day=%i', day(index)) ...
                sprintf('&time_flag=%i', ~utc) ...
                sprintf('&hour=%#g', hour(index)) ...
                '&geo_flag=' geo_flag ...
                sprintf('&latitude=%#g', latitude(index)) ...
                sprintf('&longitude=%#g', longitude(index)) ...
                sprintf('&height=%#g', altitude(index)) ...
                '&profile=1' ...
                sprintf('&start=%#g', altitude(index)) ...
                sprintf('&stop=%#g', altitude(index)) ...
                '&step=1.' ...
                endcmd]);
            if status == 0
                data((1:20) + 20*(index-1)) = parseresult(fpath, 1);
            else
                error('msis:curlError', ['Curl command did not work. ' ...
                    'It returned status:%i, ' regexprep(result, '\\', ...
                    '\\\\')], status);
            end
        end
end % End case

% data has all the data in one long vector and consists of 1 independent
% variable (height) and 19 dependent variables (the various outputs of the
% function). Therefore, a particular variable's data occurs in steps of 20
% in data. The model outputs height even though we don't care about it
% because otherwise it doesn't work in some cases.
O = data(2:20:end);
N2 = data(3:20:end);
O2 = data(4:20:end);
mass_density = data(5:20:end);
Tn = data(6:20:end);
Tex = data(7:20:end);
He = data(8:20:end);
Ar = data(9:20:end);
H = data(10:20:end);
N = data(11:20:end);
f_10_7_daily_used = data(12:20:end);
f_10_7_3month_used = data(13:20:end);
Ap_daily_used = data(14:20:end);
Ap_00_03 = data(15:20:end);
Ap_03_06 = data(16:20:end);
Ap_06_09 = data(17:20:end);
Ap_09_12 = data(18:20:end);
Ap_12_33 = data(19:20:end);
Ap_33_59 = data(20:20:end);

% Resort the output if we did a sweep.
if profile >= 1 && profile <= 7
    O = O(sortind);
    N2 = N2(sortind);
    O2 = O2(sortind);
    mass_density = mass_density(sortind);
    Tn = Tn(sortind);
    Tex = Tex(sortind);
    He = He(sortind);
    Ar = Ar(sortind);
    H = H(sortind);
    N = N(sortind);
    f_10_7_daily_used = f_10_7_daily_used(sortind);
    f_10_7_3month_used = f_10_7_3month_used(sortind);
    Ap_daily_used = Ap_daily_used(sortind);
    Ap_00_03 = Ap_00_03(sortind);
    Ap_03_06 = Ap_03_06(sortind);
    Ap_06_09 = Ap_06_09(sortind);
    Ap_09_12 = Ap_09_12(sortind);
    Ap_12_33 = Ap_12_33(sortind);
    Ap_33_59 = Ap_33_59(sortind);
end

% Output.
density = [O, N2, O2, mass_density./1e3, He, Ar, H, N] .* 1e6;
temperature = [Tn, Tex];
f_10_7_used = [f_10_7_daily_used, f_10_7_3month_used];
Ap_used = [Ap_daily_used, Ap_00_03, Ap_03_06, Ap_06_09, Ap_09_12, ...
    Ap_12_33, Ap_33_59];


% Return the day of the year.
function dayyear = dayofyear(year, month, day)

previous = cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30]);
dayyear = previous(month) + day + double( ...
    (~mod(year, 4) & mod(year, 100)) | (~mod(year, 400)) & (month > 2) );


% Parse the result file output by the system command.
function data = parseresult(fpath, sweeplen)

% Get the data from the file temp.html into a string, then delete the file.
fid = fopen(fullfile(fpath, 'temp.html'), 'r');
if fid == -1
    error('msis:parseresult:cannotOpenFile', ['Cannot open temp.html ' ...
        'file generated by curl command.']);
end
result = fread(fid, '*char').';
f = fclose(fid);
if f == -1
    warning('msis:parseresult:cannotCloseFile', ['Cannot close ' ...
        'temp.html file generated by curl command.']);
end
delete(fullfile(fpath, 'temp.html'));

% Data starts at line 37. If there isn't a line 37, there was an error.
newlines = find(result == sprintf('\n'));
if length(newlines) < 37
    % Output the error with the HTML tags removed.
    error('msis:parseresult:modelWebError', regexprep(['The online ' ...
        'interface returned:\n' result], '<[^>]*>', ''));
end
data = sscanf(result(newlines(36)+1 : newlines(36 + sweeplen)-1), '%f');