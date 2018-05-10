function amisrData = aer_to_field_aligned_coords(amisrData,projectionAlt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Input:
% projectionAlt - Altitude of projection [km]
% amisrData --> az, el, range, magBeamNo, nBeams
%
% Output:
% amisrData --> (added)
%           --> magCartCoords.xEast,yNorth,zUp % Field aligned beams
%           --> cartCoords.xEast,yNorth,zUp

% Initializing
az = double(amisrData.az);
elev = double(amisrData.el);
% slantRange = double(amisrData.range); % There is difference between range
% and altitude/sin(elev) - some numerical error perhaps in the HDF5 file
slantRange = double(amisrData.altitude./sind(elev));
magBeamNo = amisrData.magBeamNo;
nBeams = amisrData.nBeams;

% Convert AER to Cartesian
[yNorth, xEast, zDown] = aer2ned(az,elev,slantRange);
zUp = -zDown;

% Creating copies of mag. field. aligned beams
yNorthMag = repmat(yNorth(:,magBeamNo),1,nBeams);
xEastMag = repmat(xEast(:,magBeamNo),1,nBeams);
zDownMag = repmat(zDown(:,magBeamNo),1,nBeams);
zUpMag = -zDownMag;

% Interpolating Cartesian coordinates to altitude common to mag. field aligned beam
altitudeGrid = amisrData.altitude(:,magBeamNo);
hWait = waitbar(0);
for iBeam = 1:1:nBeams
    custom_waitbar(hWait,iBeam,nBeams,'Creating magnetic field aligned coordinates');
    xEast1(:,iBeam) = interp1(zUp(:,iBeam),xEast(:,iBeam),altitudeGrid,'linear','extrap');
    yNorth1(:,iBeam) = interp1(zUp(:,iBeam),yNorth(:,iBeam),altitudeGrid,'linear','extrap');
end
delete(hWait);

projectionAltNo = find_altitude(altitudeGrid,projectionAlt);
yNorthMag1 = yNorthMag - (yNorthMag(projectionAltNo)-yNorth1(projectionAltNo,:));
xEastMag1 = xEastMag - (xEastMag(projectionAltNo)-xEast1(projectionAltNo,:));

% Magnetic field aligned coordinates, interpolated to altitude points of the mag. field aligned beam
amisrData.magCartCoords.xEast = xEastMag1;
amisrData.magCartCoords.yNorth = yNorthMag1;
amisrData.magCartCoords.zUp = zUpMag;

[amisrData.magGeodeticCoords.lat,amisrData.magGeodeticCoords.lon,amisrData.magGeodeticCoords.alt]=...
    ned2geodetic(yNorthMag1,xEast,-zUpMag,amisrData.site.latitude,amisrData.site.longitude,...
    amisrData.site.altitude./1000,wgs84Ellipsoid('km'));

% Not magnetic field aligned, but interpolated to altitude points of the mag. field aligned beam
amisrData.cartCoords.xEast = xEast1;
amisrData.cartCoords.yNorth = yNorth1;
amisrData.cartCoords.zUp = zUpMag;

% Not magnetic field aligned, but interpolated to altitude points of the mag. field aligned beam
amisrData.origCartCoords.xEast = xEast;
amisrData.origCartCoords.yNorth = yNorth;
amisrData.origCartCoords.zUp = zUp;

amisrData.projectionAlt = projectionAlt;
end

