function [ra, dec] = angles (pos)

% this function converts a vector to angular quantities

% input

%  pos = position vector, equatorial rectangular coordinates

% output

%  ra  = right ascension in hours
%  dec = declination in degrees

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seccon = 206264.8062470964;

xyproj = sqrt(pos(1)^2 + pos(2)^2);

r = atan2(pos(2), pos(1));

d = atan2(pos(3), xyproj);

ra = r * seccon / 54000.0;

dec = d * seccon / 3600.0;

if (ra < 0.0)
   ra = ra + 24;
end

