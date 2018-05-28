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

% Removing altitude < 10 km
if amisrData.altitude(1,magBeamNo)<10
    minAltNo = find_altitude(amisrData.altitude(:,magBeamNo),10);
else
    minAltNo = 1;
end
maxAltNo = length(amisrData.altitude(:,magBeamNo));
altRange = minAltNo:maxAltNo; 

% Creating copies of mag. field. aligned beams
yNorthMag = repmat(yNorth(altRange,magBeamNo),1,nBeams);
xEastMag = repmat(xEast(altRange,magBeamNo),1,nBeams);
zDownMag = repmat(zDown(altRange,magBeamNo),1,nBeams);
zUpMag = -zDownMag;

% Interpolating Cartesian coordinates to altitude common to mag. field aligned beam
altitudeGrid = amisrData.altitude(altRange,magBeamNo);
projectionAltNo = find_altitude(altitudeGrid,projectionAlt);

dBeams = 1./nBeams;
multiWaitbar('Creating field-aligned coordinates',0);
for iBeam = 1:1:nBeams
    multiWaitbar('Creating field-aligned coordinates','Increment',dBeams);
    xEast1(:,iBeam) = interp1(zUp(altRange,iBeam),xEast(altRange,iBeam),altitudeGrid,'linear','extrap');
    yNorth1(:,iBeam) = interp1(zUp(altRange,iBeam),yNorth(altRange,iBeam),altitudeGrid,'linear','extrap');
end
% multiWaitbar('Creating field-aligned coordinates','Close');


yNorthMag1 = yNorthMag - (yNorthMag(projectionAltNo)-yNorth1(projectionAltNo,:));
xEastMag1 = xEastMag - (xEastMag(projectionAltNo)-xEast1(projectionAltNo,:));

% Magnetic field aligned coordinates, interpolated to altitude points of the mag. field aligned beam
amisrData.magCartCoords.xEast = xEastMag1;
amisrData.magCartCoords.yNorth = yNorthMag1;
amisrData.magCartCoords.zUp = zUpMag;

[amisrData.magGeodeticCoords.lat,amisrData.magGeodeticCoords.lon,amisrData.magGeodeticCoords.alt]=...
    ned2geodetic(amisrData.magCartCoords.yNorth,amisrData.magCartCoords.xEast,...
    -amisrData.magCartCoords.zUp,amisrData.site.latitude,amisrData.site.longitude,...
    amisrData.site.altitude./1000,wgs84Ellipsoid('km'));

% Not magnetic field aligned, but interpolated to altitude points of the mag. field aligned beam
amisrData.cartCoords.xEast = xEast1;
amisrData.cartCoords.yNorth = yNorth1;
amisrData.cartCoords.zUp = zUpMag;

% Not magnetic field aligned, but interpolated to altitude points of the mag. field aligned beam
amisrData.origCartCoords.xEast = xEast(altRange,:);
amisrData.origCartCoords.yNorth = yNorth(altRange,:);
amisrData.origCartCoords.zUp = zUp(altRange,:);

amisrData.projectionAlt = projectionAlt;

amisrData.electronDensity = amisrData.electronDensity(altRange,:,:);
amisrData.altitude = amisrData.altitude(altRange,:);
amisrData.range = amisrData.range(altRange,:);
amisrData.dNeFrac = amisrData.dNeFrac(altRange,:,:);
amisrData.az = amisrData.az(altRange,:);
amisrData.el = amisrData.el(altRange,:);


end

