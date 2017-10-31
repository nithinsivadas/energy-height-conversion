function JD = unix2julian(UT)
% converts unix time to julian dates
JD = UT/86400 +2440587.5;
