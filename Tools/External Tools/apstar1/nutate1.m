function pos2 = nutate1 (tjd, pos1)

% this function nutates equatorial rectangular coordinates from
% mean equator and equinox of epoch to true equator and equinox of
% epoch. see pages 41-45 of the explanatory supplement to the ae.

% jpl binary ephemeris and etilt1 function

% input

%  tjd  = tdb julian date of epoch
%  pos1 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean equator and equinox
%         of epoch

% output

%   pos2 = position vector, geocentric equatorial rectangular
%          coordinates, referred to true equator and equinox
%          of epoch (out)

% note: 
% if tjd is negative, inverse nutation (true to mean) is applied

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seccon = 206264.8062470964;

tjd1 = abs(tjd);

[oblm, oblt, eqeq, dpsi, deps] = etilt1 (tjd1);

oblm = oblm * 3600 / seccon;

oblt = oblt * 3600 / seccon;

dpsi = dpsi / seccon;

deps = deps / seccon;

cobm = cos(oblm);

sobm = sin(oblm);

cobt = cos(oblt);

sobt = sin(oblt);

cpsi = cos(dpsi);

spsi = sin(dpsi);

% nutation rotation matrix follows

xxn = cpsi;
yxn = -spsi * cobm;
zxn = -spsi * sobm;

xyn = spsi * cobt;
yyn = cpsi * cobm * cobt + sobm * sobt;
zyn = cpsi * sobm * cobt - cobm * sobt;

xzn = spsi * sobt;
yzn = cpsi * cobm * sobt - sobm * cobt;
zzn = cpsi * sobm * sobt + cobm * cobt;

if (tjd < 0)
   % perform inverse rotation

   pos2(1) = xxn * pos1(1) + xyn * pos1(2) + xzn * pos1(3);
   pos2(2) = yxn * pos1(1) + yyn * pos1(2) + yzn * pos1(3);
   pos2(3) = zxn * pos1(1) + zyn * pos1(2) + zzn * pos1(3);
else
   % perform rotation
      
   pos2(1) = xxn * pos1(1) + yxn * pos1(2) + zxn * pos1(3);
   pos2(2) = xyn * pos1(1) + yyn * pos1(2) + zyn * pos1(3);
   pos2(3) = xzn * pos1(1) + yzn * pos1(2) + zzn * pos1(3);
end

