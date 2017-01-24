function [ data ] = read_all_sky_fits(fileNameStr, az, el,...
    imageSize, minElevation, projectionAltitude)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data = fitsread(fileNameStr);

azOldRes = az;
elOldRes = el;

az = change_matrix_size(azOldRes, imageSize, imageSize);
el = change_matrix_size(elOldRes, imageSize, imageSize);
data = change_matrix_size(dataOldRes, imageSize, imageSize);

az_i = az(:);
el_i = el(:);
data_i = data(:);

grad_thresh = 15;
[Fx,Fy] = gradient(az);
bad_data_logic = hypot(Fx, Fy) > grad_thresh;
az(bad_data_logic) = 0;
zerodata = (az(:)==0 | el(:)<=minElevation);
keepdata = ~zerodata(:);

az_new = az_i(keepdata);
el_new = el_i(keepdata);
data_new = data_i(keepdata);

slantRange = projectionAltitude*1000./sind(el_new);

sensorloc = [65.1260,-147.4789,689 ];

coordnames = 'Spherical';
dataloc = [slantRange(:),az_new(:),el_new(:)];

[lat,lon,h] = aer2geodetic(dataloc(:,2),dataloc(:,3),dataloc(:,1),sensorloc(1),sensorloc(2),sensorloc(3),wgs84Ellipsoid('m'));


end

