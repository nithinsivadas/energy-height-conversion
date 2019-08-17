function [url] = create_DASC_url(time,wavelength)
%create_DASC_url Creates DASC url given time and wavelength
%   Detailed explanation goes here


wavelengthStr = num2str(wavelength,'%04.f');
if ~isdatetime(time)
    timeMat = unix_to_matlab_time(time);
else
    timeMat = time;
end
host = 'ftp://optics.gi.alaska.edu/PKR/DASC/RAW/';
filename = strcat('PKR_DASC_',wavelengthStr,'_',datestr(timeMat,'yyyymmdd_HHMMSS.FFF'),'.FITS');
url = string(strcat(host,datestr(timeMat,'yyyy'),'/',datestr(timeMat,'yyyymmdd'),'/',filename));
end

