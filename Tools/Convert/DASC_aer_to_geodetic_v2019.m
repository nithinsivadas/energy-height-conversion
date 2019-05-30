function [ goodDataMinElRemoved, lat, lon, alt] =...
    DASC_aer_to_geodetic_v2018(az, el,...
    minElevation, projectionAltitude, sensorLoc)
%% DASC_aer_to_geodetic Converts input Digital All Sky Camera image with az, el
% coordinates to lat, lon projected at the given projection altitude.
% This also reads from FITS file, and produce a usable matlab structure
% More accurate slant range calculation taking into account the Earth's curvature
% Need to initialize Geodata
%--------------------------------------------------------------------------
% Input:
%--------------
% fileNameStr : Name of the FITS data file
% az          : Azimuth Coordinate of data (Default form: 0 deg -> North;
%               90 deg -> East [MxN Matrix] (after calibrate_DASC_pixels
%               function)
% el          : Elevation Coordinate of data [MxN Matrix]
% minElevation: The minimum elevation angle [eg., 15]
% projectionAltitude: The atltitude where the image ought to be projected
%                     to in [km]
% sensorloc   : Specifying location of DASC
%--------------------------------------------------------------------------
% Output:
%-------------
% dataNew     : All-Sky Image data cropped within the required elevation
% lat         : Latitude points of dataNew
% lon         : Longitude points of dataNew
% az_new      : Azimuth points of dataNew, modified by removing unwanted
%               elevation
% el_new      : Elevation points of dataNew modified similarly

%
% Last Updated : 19th May 2019
% Created by Nithin Sivadas
%%
if nargin<5
    % Specifying location of the Camera
    sensorLoc = [65.1260,-147.4789,689 ]; % Pokerflat DASC
end

% %% Read FITS data from file
% dataOldRes = fitsread(fileNameStr);
% % initialize_geodata;
% %% Change resolution of data location and data to size of az & el
% data = modify_matrix_size(dataOldRes, size(az,1), size(az,2));

%% Removing Elevation < minElevation
goodDataMinElRemoved = ones(size(az));
if minElevation~=0
    zerodata = (el(:)<=minElevation);
    baddata = zerodata(:);
    goodDataMinElRemoved(baddata) = 0;
end

%% Specifying the slant range for altitude specified by ProjectAltitude
% slantRange = projectionAltitude*1000./sind(el);
RE = 6.371*10^6;
projectionAltitude1 = projectionAltitude*1000 - sensorLoc(3);
slantRange = -RE.*sind(el) + modSign(el).*sqrt((RE^2).*(sind(el)).^2 + projectionAltitude1.*(projectionAltitude1+2.*RE));
slantRange(slantRange<=projectionAltitude1*0.1) = nan;

%%
% dataloc = [slantRange(:),az_new(:),el_new(:)];
% [lat,lon,h] = aer2geodetic(dataloc(:,2),dataloc(:,3),dataloc(:,1),sensorLoc(1),sensorLoc(2),sensorLoc(3),wgs84Ellipsoid('m'));
[lat,lon,alt] = aer2geodetic(az,el,slantRange,sensorLoc(1),sensorLoc(2),sensorLoc(3),wgs84Ellipsoid('m'));

end

function y = modSign(x)
    for i=1:1:length(x)
        if x(i)==0
            x(i) = +1;
        end
    end
    y = sign(x);
end
