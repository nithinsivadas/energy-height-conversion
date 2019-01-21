% example 7: TAI-UTC

clc
clear
format long g

% TAI-UTC for 0h UTC on 2009 Feb 13.
deltat = iauDat(2009, 2, 13, 0.0);

fprintf(1,'%+5.1f\n', deltat);

