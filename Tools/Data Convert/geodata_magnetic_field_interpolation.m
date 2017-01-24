function [ pfisrGD, magcoords, atTime, atAltitude] = geodata_magnetic_field_interpolation...
    ( pfisrGD, minTimeStr, maxTimeStr, thisAltitude, magBeamNo)
%geodata_magnetic_field_interpolation Converts the data into cartesian
%coordinates, and creates beams that point along the magnetic field
%direction

% Input 
% pfisrGD      - coordinates of pfisrGD data should be in aer coordinates
% thisTimeStr  - Time string where you wish the data to be interpolated
% thisAltitude - Altitude where the beam coordinates will be 

% Output
% pfisrGD    - Data interpolated to beams parallel to mag. field line and
%              coordinates in xEast, yNorth, zUp
% magcoords  - [xEast, yNorth, zUp] 
% atTime     - Beginning and end time of the data in matlab time units
% atAltitude - The altitude where the coordinates of beams parallel to magnetic field
%              lines coincide with the actual beam coordinates

 switch nargin
        case 5
            %continue
        case 4
            magBeamNo = 1;
        case 3
            thisAltitude = 60;
            magBeamNo = 1;
        case 2
            maxTimeStr = '26 Mar 2008 11:48:56';
            thisAltitude = 60; % km
            magBeamNo = 1;
        case 1
            maxTimeStr = '26 Mar 2008 11:48:56';
            minTimeStr = '26 Mar 2008 11:48:56';
            thisAltitude = 60; % km
            magBeamNo = 1;
        otherwise
            if nargin>4
                error('TooManyInputs: Number of input variables should be at most 4')
            else
                error('TooFewInputs: Number of input variables should at least be 1');
            end
    end

% Converting aer (azimuth, elevation, slant height) to Cartesian (North
% East Down)
dataLocation = pfisrGD.dataloc;
az   = dataLocation(:,2);
elev = dataLocation(:,3);
slant= dataLocation(:,1);
[yNorth, xEast, zDown] = aer2ned(az,elev,slant);
 zUp = -zDown;
 
 % In JSwoboda's geodata class, the cartesian coordinates use East North Up
 % coordinates. That is, the xEast, yNorth, zUp. 

% Creating copies of the magnetically field aligned beams 
azMag   = repmat(az(magBeamNo),length(az),1);
elevMag = repmat(elev(magBeamNo ),length(elev),1);
slantMag= slant;

[yMagNorth, xMagEast, zMagDown] = aer2ned(azMag,elevMag,slantMag);
 zMagUp = -zMagDown;

 % Calculating the number of beams
coordinateNo=1; % Using slant range to evaluate the number of beams
nBeams = calculate_nBeams(pfisrGD, coordinateNo);
 
% Creating the default altitude grid (altitude grid of magnetic field
% aligned beam)
altitudeGrid = zMagUp(magBeamNo:nBeams:end);
altitudeNo = find_altitude( altitudeGrid, thisAltitude);
atAltitude = altitudeGrid(altitudeNo);
 
% Interpolating the cartesian coordinates of data locations to altitude
% common to the magnetic field aligned beam
 for i=1:1:nBeams
     xEast(i:nBeams:end) = interp1(zUp(i:nBeams:end),xEast(i:nBeams:end),altitudeGrid,'linear','extrap');
     yNorth(i:nBeams:end) = interp1(zUp(i:nBeams:end),yNorth(i:nBeams:end),altitudeGrid,'linear','extrap');
 end;
 
 % Forming coordinates of beams that are aligned with the magnetic field line
 for i=1:1:length(yMagNorth)
     BeamNo=rem(i,nBeams)+altitudeNo*nBeams; % Corresponding to an altitude of 60.55 km
     if rem(i,nBeams)==0
         BeamNo=nBeams+altitudeNo*nBeams;
     end;
     y1MagNorth(i) = yMagNorth(i)-(yMagNorth(BeamNo)-yNorth(BeamNo));
     x1MagEast(i)  = xMagEast(i)-(xMagEast(BeamNo)-xEast(BeamNo)); 
 end;
 
z1MagUp  = zMagUp;

minTimeNo = find_time(unix_to_matlab_time(pfisrGD.times(:,1)),minTimeStr);
maxTimeNo = find_time(unix_to_matlab_time(pfisrGD.times(:,1)),maxTimeStr);
if minTimeNo==maxTimeNo
    timeNo = minTimeNo;
elseif minTimeNo>maxTimeNo
    error('Lower limit of time is greater than the upper limit');
else
    timeNo = minTimeNo:1:maxTimeNo;
end;
atTime = unix_to_matlab_time(pfisrGD.times(timeNo,:));
pfisrGD.timereduce(timeNo);

% Interpolating the data in pfisrGD in cartesian coordinates, of beams
% aligned along magnetic field
magcoords = [x1MagEast(:), y1MagNorth(:), z1MagUp(:)];

pfisrGD.interpolate(magcoords,'Cartesian','natural');

end

