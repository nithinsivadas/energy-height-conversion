function pos2 = precess (tjd1, pos1, tjd2)

% this function precesses equatorial rectangular coordinates from
% one epoch to another.  the coordinates are referred to the mean
% equator and equinox of the two respective epochs. see pages 30-34
% of the explanatory supplement to the ae, lieske, et al. (1977)
% astronomy and astrophysics 58, 1-16, and lieske (1979) astronomy
% and astrophysics 73, 282-284.

% input

%  tjd1 = tdb julian date of first epoch

%  pos1 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean equator and equinox of
%         first epoch
%  tjd2 = tdb julian date of second epoch

% output

%  pos2 = position vector, geocentric equatorial rectangular
%         coordinates, referred to mean equator and equinox of
%         second epoch

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global tlastp1 tlastp2
 
% global xxp xyp xzp yxp yyp yzp zxp zyp zzp

seccon = 206264.8062470964;

% if (abs(tjd1 - tlastp1) < 1e-6 && abs(tjd2 - tlastp2) < 1e-6)
%    pos2(1) = xxp * pos1(1) + yxp * pos1(2) + zxp * pos1(3);
%    pos2(2) = xyp * pos1(1) + yyp * pos1(2) + zyp * pos1(3);
%    pos2(3) = xzp * pos1(1) + yzp * pos1(2) + zzp * pos1(3);
%    
%    return;
% end
 
% if (abs(tjd1 - tlastp2) < 1e-6 && abs(tjd2 - tlastp1) < 1e-6)
%    % perform inverse rotation
% 
%    pos2(1) = xxp * pos1(1) + xyp * pos1(2) + xzp * pos1(3);
%    pos2(2) = yxp * pos1(1) + yyp * pos1(2) + yzp * pos1(3);
%    pos2(3) = zxp * pos1(1) + zyp * pos1(2) + zzp * pos1(3);
%    
%    return;
% end

% t0 and t below correspond to lieske's big t and little t

t0 = (tjd1 - 2451545) / 36525;

t = (tjd2 - tjd1) / 36525;

t02 = t0 * t0;

t2 = t * t;

t3 = t2 * t;

% zeta0, zee, and theta below correspond to lieske's zeta-sub-a,
% z-sub-a, and theta-sub-a

zeta0 = (2306.2181 + 1.39656 * t0 - 0.000139*t02) * t ...
        + (0.30188 - 0.000344 * t0) * t2 +  0.017998 * t3;

zee = (2306.2181 + 1.39656 * t0 - 0.000139 * t02) * t ...
      + (1.09468 + 0.000066 * t0) * t2 + 0.018203 * t3;

theta = (2004.3109 - 0.85330 * t0 - 0.000217 * t02) * t ...
        + (-0.42665 - 0.000217 * t0) * t2 - 0.041833 * t3;

zeta0 = zeta0 / seccon;

zee = zee / seccon;

theta = theta / seccon;

czeta0 = cos(zeta0);

szeta0 = sin(zeta0);

czee = cos(zee);

szee = sin(zee);

ctheta = cos(theta);

stheta = sin(theta);

% precession rotation matrix follows

xxp = czeta0 * ctheta * czee - szeta0 * szee;
yxp = -szeta0 * ctheta * czee - czeta0 * szee;
zxp = -stheta * czee;

xyp = czeta0 * ctheta * szee + szeta0 * czee;
yyp = -szeta0 * ctheta * szee + czeta0 * czee;
zyp = -stheta * szee;

xzp = czeta0 * stheta;
yzp = -szeta0 * stheta;
zzp = ctheta;

% tlastp1 = tjd1;
 
% tlastp2 = tjd2;

% perform rotation

pos2(1) = xxp * pos1(1) + yxp * pos1(2) + zxp * pos1(3);
pos2(2) = xyp * pos1(1) + yyp * pos1(2) + zyp * pos1(3);
pos2(3) = xzp * pos1(1) + yzp * pos1(2) + zzp * pos1(3);
