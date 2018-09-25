function [flag,probeNames] = find_geodetic_conjunction(probes,conjunctionCoords,stopAlt)
%% find_geodetic_conjunction.m Calculates conjunction between geodetic coordinates at stopAlt heights
%   Step 1: Calculate geodetic distance
%   Step 2: If distance is within the specified radius then flag = 1 else 0
%   --------------------------------------------
%   Input
%   ------
%   probes.(field).GDZ       - [lat, lon, alt] of a spacecraft
%   conjunctionCoords.GDZ    - [lat, lon, alt] of probe/ground-location with which you wish to check conjunction.
%   conjunctionCoords.radius - [in km] The acceptable radius of conjunction
%   stopAlt                  - [in km] Altitude at which conjunction is to be determined
%   --------------------------------------------
%   Output
%   ------
%   flag    "1" - conjunction
%           "0" - no conjunction
%   probeNames  - fieldnames of probe.(field) 
%   ----------------------------------------------

%% Defining RE
RE = (6371+stopAlt)*1000; %[m]
probeNames = fieldnames(probes);
nProbes = length(probeNames);

for thisProbe = 1:1:nProbes
    distance=geodetic_distance(probes.(char(probeNames(thisProbe))).GDZ,...
        conjunctionCoords.GDZ, RE);
    decision(:,thisProbe)=(distance<conjunctionCoords.radius.*1000);
end
flag = decision;

probeNames = probeNames';
end

