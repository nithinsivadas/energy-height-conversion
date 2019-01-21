% example 5:
% Here is example code that takes a time expressed as TAI, encodes it into
% the internal format using the DTF2D routine, and transforms it into UTC;
% the TAI and UTC are each decoded using the D2DTF routine and reported
% Encode TAI date and time into internal format.

clc
clear
format long g

[a1, a2] = iauDtf2d('TAI', 2009, 1, 1, 0, 0, 33.7);

% Decode and report the TAI.
[iy, im, id, ihmsf] = iauD2dtf('TAI', 3, a1, a2);

fprintf(1,'TAI %4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d\n',iy,im,id,ihmsf(1),...
        ihmsf(2),ihmsf(3),ihmsf(4));

% Transform TAI to UTC.
[u1, u2] = iauTaiutc(a1, a2);

% Decode and report the UTC.
[iy, im, id, ihmsf] = iauD2dtf('UTC', 3, u1, u2);

fprintf(1,'UTC %4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d\n',iy,im,id,ihmsf(1),...
        ihmsf(2),ihmsf(3),ihmsf(4));
fprintf(1,'\n');

