% example 6: transform UTC into other times

clc
clear
format long g

% Site terrestrial coordinates (WGS84).
latnd = 19;
latnm = 28;
slatn = 52.5;
lonwd = 155;
lonwm = 55;
slonw = 59.6;
hm = 0.0;

% Transform to geocentric.
phi = iauAf2a('+', latnd, latnm, slatn);
elon = iauAf2a('-', lonwd, lonwm, slonw);
xyz = iauGd2gc(1, elon, phi, hm);
u = sqrt ( xyz(1)*xyz(1) + xyz(2)*xyz(2) );
v = xyz(3);

% UTC date and time.
iy = 2006;
mo = 1;
id = 15;
ih = 21;
im = 24;
sec = 37.5;

% Transform into internal format.
[utc1, utc2] = iauDtf2d('UTC', iy, mo, id, ih, im, sec);

% UT1-UTC (s, from IERS).
dut = 0.3341;

% UTC -> UT1.
[ut11, ut12] = iauUtcut1(utc1, utc2, dut);

% Extract fraction for TDB-TT calculation, later.
ut = mod( mod(ut11,1.0) + mod(ut12,1.0), 1.0 );

% UTC -> TAI -> TT -> TCG.
[tai1, tai2] = iauUtctai(utc1, utc2);
[tt1, tt2] = iauTaitt(tai1, tai2);
[tcg1, tcg2] = iauTttcg(tt1, tt2);

% TDB-TT (using TT as a substitute for TDB).
dtr = iauDtdb ( tt1, tt2, ut, elon, u, v );
% TT -> TDB -> TCB.
[tdb1, tdb2] = iauTttdb(tt1, tt2, dtr);

[tcb1, tcb2] = iauTdbtcb(tdb1, tdb2);

% Report.
[iy, mo, id, ihmsf] = iauD2dtf('UTC', 6, utc1, utc2);

fprintf(1,'UTC%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('ut1', 6, ut11, ut12);

fprintf(1,'UT1%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('tai', 6, tai1, tai2);

fprintf(1,'TAI%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('tt', 6, tt1, tt2);

fprintf(1,'TT %5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('tcg', 6, tcg1, tcg2);

fprintf(1,'TCG%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('tdb', 6, tdb1, tdb2);

fprintf(1,'TDB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
[iy, mo, id, ihmsf] = iauD2dtf('tcb', 6, tcb1, tcb2);

fprintf(1,'TCB%5d/%2.2d/%2.2d%3d:%2.2d:%2.2d.%6.6d\n',...
iy, mo, id , ihmsf(1), ihmsf(2), ihmsf(3), ihmsf(4) );
fprintf(1,'\n');

