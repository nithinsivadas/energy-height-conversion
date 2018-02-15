function [flag,probeNames] = find_geodetic_conjunction(probes,conjunctionCoords,stopAlt)
%% find_geodetic_conjunctoin.m Calculates conjunction between geodetic coordinates at stopAlt heights
%   Detailed explanation goes here
%   probes.(field).GDZ (geodetic coordinates)
%   conjunctionCoords.GDZ
%   conjunctionCoords.radius - in km
%   stopAlt - in km (altitude at which conjunction is to be determined)
%   can add more decision parameters

%%
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

