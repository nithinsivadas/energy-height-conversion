function [az,el] = get_star_az_el...
    (RA,DEC,pmRA,pmDEC,parallax,RV,time,glat,glong,galt,...
    temperature,pressure,humidity,wavelength)
%get_star_az_el Get stars azimuth and elevation in the sky given the input star catalogue values.
%   Detailed explanation goes here
% RA,DA         double  ICRS right ascension & declination at J1991.25 (radians, Note 1)
% pmRA          double  RA proper motion (radians/year; Note 2)
% pmDEC         double  Dec proper motion (radians/year)
% parallax      double  parallax (arcsec)
% RV            double  radial  velocity (km/s, positive if receding)
% glong         double  longitude (radians, east-positive, Note 6)
% glat          double  geodetic latitude (radians, Note 6)
% galt          double  height above ellipsoid (m, geodetic Notes 6,8)
% temperature   double  ambient temperature at the observer (deg K)
% pressure      double  pressure at the observer (hPa = mB, Note 8)
% humidity      double  relative humidity at the observer (range 0-1)
% wavelength    double  wavelength (micrometers, Note 9)

%% To be done
% Extract polar motion coordinates from IERS
%
%%
t = datevec(time);
[jd_utc1,jd_utc2] = iauDtf2d('UTC', t(1), t(2), t(3), t(4), t(5), t(6));

% height above ellipsoid? what is the difference from altitude?
if nargin<11 || isempty(temperature)
   t0 = 256; %[k] - Alaska temperature
else
   t0 = temperature; %ambient temperature
end
if nargin<12 || isempty(pressure)
    pressure = 1013.25*exp(-galt./29.3*t0); % hPa
end

tempC = t0-273;

if nargin<13 || isempty(humidity)
    humidity = 0.5;
end

if nargin<14 || isempty(wavelength)
    wavelength = 557.7*10^-3; % nm to micrometers
end

dut1 = 0; % for now, get from table IERS bulletin
xp = 0; %polar motion coordinates in radians
yp = 0; %polar motion coordinates in radians

arrL = length(RA);

multiWaitbar('Creating Sky Chart...',0);
id = 1./arrL;
for i = 1:1:arrL
multiWaitbar('Creating Sky Chart...','Increment',id);
[az(i), zenithDistance(i), hourAngle(i), declination(i), rightAscension(i)]...
    = iauAtco13(RA(i),DEC(i),pmRA(i),pmDEC(i),...
    parallax(i), RV(i), jd_utc1, jd_utc2,...
    dut1, glong, glat,...
    galt, xp, yp,...
    pressure, tempC, ...
    humidity, wavelength);
end
multiWaitbar('Close All');
az = rad2deg(az);
el = 90-rad2deg(zenithDistance);

end
