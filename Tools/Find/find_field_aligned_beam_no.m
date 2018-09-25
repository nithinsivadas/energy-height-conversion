function [magBeamNo] = find_field_aligned_beam_no(beamCodes,magFieldAlignedAz,magFieldAlignedElev)
%find_field_aligned_beam_no.m Finds the field aligned beam number given the
%field aligned Az and Elev values, using the beamCodes matrix from
%SRI-AMISR data sets
%------------------------------------------------------------------
% Input
%--------
% beamCodes             - [nPos x nBeams] Beam code for each azimuth and elevation coordinate.
% magFieldAlignedAz     - The azimuth value of the field line from the radar center
% magFieldAlignedElev   - The elevation value of the field line from the radar centre
%-------------------------------------------------------------------
% Output
%--------
% magBeamNo             - Field aligned beam number
%-------------------------------------------------------------------
% 24 Jan 2018

%   Detailed explanation goes here

indxElev = find(abs(beamCodes(3,:)-magFieldAlignedElev)<=min(abs(beamCodes(3,:)-magFieldAlignedElev)));
indxAz = find(abs(beamCodes(2,:)-magFieldAlignedAz)<=min(abs(beamCodes(2,:)-magFieldAlignedAz)));

combinedIndx = [indxElev,indxAz];
[magBeamNo, times] = mode(combinedIndx);

if times<2
    error('Could not find the magnetic field aligned beam number since the number of time and Index repeates is less than 2.');
end

end

