function [ dataNew, lat, lon, az_new, el_new, sensorloc, time ] = DASC_aer_to_geodetic(fileNameStr, az, el,...
    imageSize, minElevation, projectionAltitude)
%% DASC_aer_to_geodetic Converts input Digital All Sky Camera image with az, el coordinates to lat, lon projected at the given projection altitude. 
% Input:
% fileNameStr : Name of the FITS data file
% az          : Azimuth Coordinate of data (Default form: 0 deg -> West;
%               90 deg -> North [MxN Matrix]
% el          : Elevation Coordinate of data [MxN Matrix]
% imageSize   : The size to which you would like the FITS data to be
%               interpolated to [Eg., 512]
% minElevation: The minimum elevation angle [eg., 15]
% projectionAltitude: The atltitude where the image ought to be projected
%                     to in [km]

% Output:
% dataNew     : All Sky Image data cropped within the required elevation [1-D Array]
% lat         : Latitude points of dataNew [1-D Array] 
% lon         : Longitude points of dataNew [1-D Array]
% az_new      : Azimuth points of dataNew [1-D Array]; 0 deg -> West 
%               and 90 -> North
% el_new      : Elevation points of dataNew [1-D Array] 
% sensorloc   : Specifying location of DASC

% Last Updated : 7th Feb 2017
% Created by Nithin Sivadas
%% Read FITS data from file
dataOldRes = fitsread(fileNameStr);
azOldRes = az;
elOldRes = el;

%% Change resolution of data location and data to user specified imageSize
az = modify_matrix_size(azOldRes, imageSize, imageSize);
el = modify_matrix_size(elOldRes, imageSize, imageSize);
data = modify_matrix_size(dataOldRes, imageSize, imageSize);

%% Fix problems with the coordinate matrix
% Look for large gradients in the az mapping because the in between values
% will put the data in the wrong spot. [From John Swoboda]

grad_thresh = 15;
[Fx,Fy] = gradient(az);
bad_data_logic = hypot(Fx, Fy) > grad_thresh;
az(bad_data_logic) = 0;

%% Removing Zeros and Elevation < minElevation
zerodata = (az(:)==0 | el(:)<=minElevation);
keepdata = ~zerodata(:);

az_new = az(keepdata);
el_new = el(keepdata);
dataNew = data(keepdata);

%% Specifying the slant range for altitude specified by ProjectAltitude
slantRange = projectionAltitude*1000./sind(el_new);

%% Specifying location of the Camera
sensorloc = [65.1260,-147.4789,689 ];

%% Generating Time Stamp
aldtnum = fitsfiletimestamp(fileNameStr);
time = (aldtnum-datenum('jan-01-1970'))*(24*3600);
time = unix_to_matlab_time(time);


%% 
coordnames = 'Geodetic';
dataloc = [slantRange(:),rotate_array(az_new(:),-90),el_new(:)];

% Rotating azimuth counter-clockwise by 90 degrees, as calibration data is stored
% 90 deg clockwise: East -> 0 deg, North ->270 deg
% Rotating => North -> 0 deg, and East -> 90 deg
[lat,lon,h] = aer2geodetic(dataloc(:,2),dataloc(:,3),dataloc(:,1),sensorloc(1),sensorloc(2),sensorloc(3),wgs84Ellipsoid('m'));


end

