function [data] = read_amisr(fileNameStr)
% read_amisr.m Reads the specific data in h5 format obtained from amisr.com
% -------------------------------------------------------------------------
% Input
% -------------------------------------------------------------------------
% fileNameStr : File path and name of the amisr h5 file
% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
% data.electronDensity  : Electron density profile with no Tr [nhxnT] [m-3] 
% data.altitude         : Altitude array [nhx1] [km]
% data.range            : Range array [km]
% data.time             : Time array [nTx1] [matlab units]
% data.dNeFrac          : Fractional error in Ne [nhxnT]
% data.site.latitude    : Latitude of the instrument [deg]
% data.site.longitude   : Longitude of the instrument [deg]
% data.site.altitude    : Altitude of the instrument [m]
% data.magBeamNo        : If instrument is PFISR, identifies the magnetic field aligned beam number
% data.nBeams           : Total number of beams available
% data.az               : Azimuth of the data points
% data.el               : Elevation of hte data points
%---------------------------------------------------------------
% Created by: Nithin Sivadas
% Date      : 
%---------------------------------------------------------------

data.electronDensity = h5read(fileNameStr,'/NeFromPower/Ne_NoTr'); % in [m^-3]
data.altitude = h5read(fileNameStr,'/NeFromPower/Altitude')/1000; % in [km]
data.range = h5read(fileNameStr,'/NeFromPower/Range')/1000; % in [km]
data.time = unix_to_matlab_time(double(h5read(fileNameStr,'/Time/UnixTime'))); % [matlab time]
% data.electronDensity(data.electronDensity<0)=10^6; %this isn't a good idea
data.dNeFrac = h5read(fileNameStr,'/NeFromPower/dNeFrac'); % Fractional error in Ne
data.site.latitude = h5read(fileNameStr,'/Site/Latitude');
data.site.longitude = h5read(fileNameStr,'/Site/Longitude');
data.site.altitude = h5read(fileNameStr,'/Site/Altitude'); % [in m]
data.site.code = h5read(fileNameStr,'/Site/Code');
beamCodes = h5read(fileNameStr,'/BeamCodes');

if data.site.code == 61
    data.magBeamNo = find_field_aligned_beam_no(beamCodes,-154.3,77.5);
end

nBeams = (size(beamCodes,2));
for iBeam=1:1:nBeams
    data.az(:,iBeam) = ones(size(data.altitude(:,iBeam))).*(beamCodes(2,iBeam));
    data.el(:,iBeam) = ones(size(data.altitude(:,iBeam))).*(beamCodes(3,iBeam));
end
% data.slantHeight = data.altitude./sind(data.el);
data.range = repmat(data.range,1,nBeams);
data.nBeams = nBeams;

end

