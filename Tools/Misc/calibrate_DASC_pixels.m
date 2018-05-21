function [az, el, goodData] = calibrate_DASC_pixels(azCal,elCal,imageSize,dateCal)
%calibrate_DASC_pixels Use data from calibration FIT files to calculate the
%az, el of each pixel - as well as cut off bad data with large az
%gradients.
%--------------------------------------------------------------------------
% Input:
%--------------
% azCal       : Azimuth Coordinate of data from cal file (Default form: 0 deg -> West;
%               90 deg -> North [MxN Matrix]
% elCal       : Elevation Coordinate of data from cal file [MxN Matrix]
% imageSize   : The size to which you would like the FITS data to be
%               interpolated to [Eg., 512]
% dateCal     : The matlab datenum value of the calibration file
%--------------------------------------------------------------------------
% Output:
%-------------
% az          : Azimuth points rotated in the right direction (0 deg -> North 
%               and 90 -> East)
% el          : Elevation points after rotating [1-D Array] 
% goodData    : Logic array, with all the good coordinates labeled 1, and
%               bad coordinates (with an grad(az)>15) labelled 0.
% Last Updated : 18th May 2018
% Created by Nithin Sivadas

if nargin<4
    dateCal = datenum('20111006','yyyymmdd');
end
if nargin<3
    % Specifying location of the Camera
    imageSize = 512;
end

%% Change resolution of data location and data to user specified imageSize
az = modify_matrix_size(azCal, imageSize, imageSize);
el = modify_matrix_size(elCal, imageSize, imageSize);

%% Fix problems with the coordinate matrix
% Look for large gradients in the az mapping because the in between values
% will put the data in the wrong spot. [From John Swoboda]
grad_thresh = 15;
[Fx,Fy] = gradient(az);
bad_data_logic = hypot(Fx, Fy) > grad_thresh;
goodData=ones(size(az));
goodData(bad_data_logic) = 0;


%% Rotate az, and el of cal file to appropriate orientation
if dateCal<=datenum('20111006','yyyymmdd')
    az = rotate_array(az,-90);
    % Rotating azimuth counter-clockwise by 90 degrees, as calibration data is stored
    % 90 deg clockwise: East -> 0 deg, North ->270 deg
    % Rotating => North -> 0 deg, and East -> 90 deg
else
    error('Do not know how to rotate the image for these data ranges');
end

end

