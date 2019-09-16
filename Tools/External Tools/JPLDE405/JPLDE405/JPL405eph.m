clc
clear all
format long g

SAT_Const

global PC

DE405

% Initialize UT1-UTC and TAI-UTC time difference
fid = fopen('eop19620101.txt','r');

%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------

eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);

fclose(fid);

Mjd_UTC = Mjday(2008,8,1,0,0,0);

[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

Mjd_TT = Mjd_UTC + TT_UTC/86400.0;

[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE405(Mjd_TT)

