% example 1: UTC to TT

clc
clear
format long g

% Encode UTC date and time into internal format.
[u1, u2] = iauDtf2d ('UTC', 2010, 7, 24, 11, 18, 7.318);

% Transform UTC to TAI, then TAI to TT.
[a1, a2] = iauUtctai(u1, u2);
[t1, t2] = iauTaitt(a1, a2);

% Decode and report the TT.
[iy, im, id, ihmsf] = iauD2dtf('tt', 3, t1, t2);

fprintf(1,'%4d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%3.3d\n',iy,im,id,ihmsf(1),ihmsf(2),ihmsf(3),ihmsf(4));
fprintf(1,'\n');

