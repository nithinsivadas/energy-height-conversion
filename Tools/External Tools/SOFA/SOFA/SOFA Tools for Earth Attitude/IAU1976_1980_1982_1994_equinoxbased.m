% IAU 2006/2000A, CIO based, using classical angles

clc
clear all
format long g

constants

% UTC
IY = 2007;
IM = 4;
ID = 5;
IH = 12;
MIN = 0;
SEC = 0;

% Polar motion (arcsec->radians).
XP = 0.0349282 * AS2R;
YP = 0.4833163 * AS2R;
% UT1-UTC (s).
DUT1 = -0.072073685;
% Nutation corrections wrt IAU 1976/1980 (mas->radians).
DDP80 = -55.0655 * AS2R/1000;
DDE80 = -6.3580 * AS2R/1000;
% % CIP offsets wrt IAU 2000A (mas->radians).
% DX00 = 0.1725 * AS2R/1000;
% DY00 = -0.2650 * AS2R/1000;
% % CIP offsets wrt IAU 2006/2000A (mas->radians).
% DX06 = 0.1750 * AS2R/1000;
% DY06 = -0.2259 * AS2R/1000;

% TT (MJD)
[DJMJD0 DATE] = iauCal2jd(IY, IM, ID);
TIME = ( 60*(60*IH + MIN) + SEC ) / 86400;
UTC = DATE + TIME;
DAT = iauDat ( IY, IM, ID, TIME );
TAI = UTC + DAT/86400;
TT = TAI + 32.184/86400;
% UT1
TUT = TIME + DUT1/86400;
UT1 = DATE + TUT;
% IAU 1976 precession matrix, J2000.0 to date.
RP = iauPmat76(DJMJD0, TT);
% IAU 1980 nutation.
[DP80, DE80] = iauNut80(DJMJD0, TT);
% Add adjustments: frame bias, precession-rates, geophysical.
DPSI = DP80 + DDP80;
DEPS = DE80 + DDE80;
% Mean obliquity.
EPSA = iauObl80(DJMJD0, TT);
% Build the rotation matrix.
RN = iauNumat(EPSA, DPSI, DEPS);
% Combine the matrices: PN = N x P.
RNPB = RN * RP
% Equation of the equinoxes, including nutation correction.
EE = iauEqeq94 ( DJMJD0, TT ) + DDP80 * cos ( EPSA );
% Greenwich apparent sidereal time (IAU 1982/1994).
GST = iauAnp ( iauGmst82 ( DJMJD0+DATE, TUT ) + EE )
% Form celestial-terrestrial matrix (no polar motion yet).
RC2TI = iauCr(RNPB);
RC2TI = iauRz(GST, RC2TI)
% Polar motion matrix (TIRS->ITRS, IERS 1996).
RPOM = zeros(3);
RPOM = iauIr(RPOM);
RPOM = iauRx(-YP, RPOM);
RPOM = iauRy(-XP, RPOM);
% Form celestial-terrestrial matrix (including polar motion).
RC2IT = iauRxr(RPOM, RC2TI)

