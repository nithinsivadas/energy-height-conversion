% Example:
% Starting with catalog data for a star, the demonstration program listed 
% below performs a whole series of transformations, namely:
% 1. ICRS to CIRS.
% 2. The reverse, giving the astrometric place.
% 3. Astrometric place to CIRS.
% 4. Geocentric apparent place via the equation of the origins.
% 5. CIRS to topocentric, i.e. CIRS to observed but with zero air pressure.
% 6. CIRS to observed.
% 7. ICRS to observed in a single call.
% 8. ICRS to CIRS using JPL DE405 for the Earth ephemeris.
% 9. The same but including light de ection by Jupiter and Saturn.
% 10. The reverse, to check agreement with Step 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long g

constants

% Site longitude, latitude (radians) and height above the geoid (m).
elong = iauAf2a ('-', 5, 41, 54.2);
phi =  iauAf2a ('-', 15, 57, 42.8);
hm = 625.0;
% Ambient pressure (HPa), temperature (C) and rel. humidity (frac).
phpa = 952.0;
tc = 18.5;
rh = 0.83;
% Effective color (microns).
wl = 0.55;
% UTC date.
[utc1, utc2] = iauDtf2d('UTC', 2013, 4, 2, 23, 15, 43.55);
% TT date.
[tai1, tai2] = iauUtctai (utc1, utc2);
[tt1, tt2] = iauTaitt (tai1, tai2);
% EOPs: polar motion in radians, UT1-UTC in seconds.
xp = 50.995e-3 * DAS2R;
yp = 376.723e-3 * DAS2R;
dut1 = 155.0675e-3;
% Corrections to IAU 2000A CIP (radians).
dx = 0.269e-3 * DAS2R;
dy = -0.274e-3 * DAS2R;
% Star ICRS RA,Dec (radians).
rc = iauTf2a (' ', 14, 34, 16.81183);
dc = iauAf2a ('-', 12, 31, 10.3965);
fprintf (1,'     ICRS, epoch J2000.0:');
[~, i] = iauA2tf(7, rc);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, dc);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% Proper motion: RA/Dec derivatives, epoch J2000.0.
pr = atan2 (-354.45e-3 * DAS2R, cos(dc));
pd = 595.35e-3 * DAS2R;
% Parallax (arcsec) and recession speed (km/s).
px = 164.99e-3;
rv = 0.0;
% ICRS to CIRS (geocentric observer).
[ri, di, eo] = iauAtci13 (rc, dc, pr, pd, px, rv, tt1, tt2);
fprintf (1,'         catalog -> CIRS:');
[~, i] = iauA2tf(7, ri);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, di);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% CIRS to ICRS (astrometric).
[rca, dca, eo] = iauAtic13 (ri, di, tt1, tt2);
fprintf (1,'     CIRS -> astrometric:');
[~, i] = iauA2tf(7, rca);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, dca);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% ICRS (astrometric) to CIRS (geocentric observer).
[ri, di, eo] = iauAtci13 (rca, dca, 0.0, 0.0, 0.0, 0.0, tt1, tt2);
fprintf (1,'     astrometric -> CIRS:');
[~, i] = iauA2tf(7, ri);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, di);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% Apparent place.
ra = iauAnp (ri - eo);
da = di;
fprintf (1,'     geocentric apparent:');
[~, i] = iauA2tf(7, ra);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, da);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% CIRS to topocentric.
[aot, zot, hot, dot, rot] = iauAtio13(ri, di, utc1, utc2, dut1, elong, ...
                                      phi, hm, xp, yp, 0.0, 0.0, 0.0, 0.0);
fprintf (1,'     CIRS -> topocentric:');
[~, i] = iauA2tf(7, rot);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, dot);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% CIRS to observed.
[aob, zob, hob, dob, rob] = iauAtio13(ri, di, utc1, utc2, dut1, elong, ...
                                      phi, hm, xp, yp, phpa, tc, rh, wl);
fprintf (1,'        CIRS -> observed:');
[~, i] = iauA2tf(7, rob);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, dob);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% ICRS to observed.
[aob, zob, hob, dob, rob, eo] = iauAtco13(rc, dc, pr, pd, px, rv, utc1,...
                     utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl);
fprintf (1,'        ICRS -> observed:');
[~, i] = iauA2tf(7, rob);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af(6, dob);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% ICRS to CIRS using some user-supplied parameters.
% SOFA heliocentric Earth ephemeris.
[pvh, pvb] = iauEpv00 (tt1, tt2);
% JPL DE405 barycentric Earth ephemeris.
pvb(1,1) = -0.9741704366519668;
pvb(1,2) = -0.2115201000882231;
pvb(1,3) = -0.0917583114068277;
pvb(2,1) = 0.0036436589347388;
pvb(2,2) = -0.0154287318503146;
pvb(2,3) = -0.0066892203821059;
% IAU 2000A CIP.
r = iauPnm00a (tt1, tt2);
[x, y] = iauBpn2xy (r);
% Apply IERS corrections.
x = x + dx;
y = y + dy;
% SOFA CIO locator.
s = iauS06 ( tt1, tt2, x, y );
% Populate the context.
astrom = iauApci (tt1, tt2, pvb, pvh(1,:), x, y, s);
% Carry out the transformation and report the results.
[ri, di] = iauAtciq ( rc, dc, pr, pd, px, rv, astrom);
fprintf (1,'ICRS -> CIRS (JPL, IERS):');
[~, i] = iauA2tf(7, ri);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, di);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% The same but with Saturn then Jupiter then Sun light deflection.
b(1).bm = 0.00028574;
b(1).dl = 3e-10;
b(1).pv(1,1) = -7.8101442680818964;
b(1).pv(1,2) = -5.6095668114887358;
b(1).pv(1,3) = -1.9807981923749924;
b(1).pv(2,1) = 0.0030723248971152;
b(1).pv(2,2) = -0.0040699547707598;
b(1).pv(2,3) = -0.0018133584165345;
b(2).bm = 0.00095435;
b(2).dl = 3e-9;
b(2).pv(1,1) = 0.7380987962351833;
b(2).pv(1,2) = 4.6365869247538951;
b(2).pv(1,3) = 1.9693136030111202;
b(2).pv(2,1) = -0.0075581692172088;
b(2).pv(2,2) = 0.0012691372216750;
b(2).pv(2,3) = 0.0007279990012801;
b(3).bm = 1.0;
b(3).dl = 6e-6;
b(3).pv(1,1) = -0.0007121743770509;
b(3).pv(1,2) = -0.0023047830339257;
b(3).pv(1,3) = -0.0010586596574639;
b(3).pv(2,1) = 0.0000062923521264;
b(3).pv(2,2) = -0.0000003308883872;
b(3).pv(2,3) = -0.0000002964866231;
[ri, di] = iauAtciqn(rc, dc, pr, pd, px, rv, astrom, 3, b);
fprintf (1,'ICRS -> CIRS (+ planets):');
[~, i] = iauA2tf(7, ri);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, di);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

% CIRS to ICRS (astrometric).
[rca, dca] = iauAticqn(ri, di, astrom, 3, b);
fprintf (1,'     CIRS -> astrometric:');
[~, i] = iauA2tf(7, rca);
fprintf(1,' %2.2d %2.2d %2.2d.%7.7d', i(1),i(2),i(3),i(4));
[pm, i] = iauA2af (6, dca);
fprintf(1,' %c%2.2d %2.2d %2.2d.%6.6d\n', pm,i(1),i(2),i(3),i(4));

