% demo_astar1.m    April 28, 2008

% this script demonstrates how to use the apstar1
% Matlab function to calculate the apparent geocentric
% or topocentric coordinates of a star

% jpl binary ephemeris

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global iephem km ephname

rtd = 180 / pi;

% name of binary ephemeris file

ephname = 'de421.bin';

% initialize jplephem

iephem = 1;

% request output in au and au/day

km = 0;

% read leap seconds data file

readleap;

% conversion factor - astronomical unit to kilometers

au = 149597870.691;

% request calendar date

clc; home;

fprintf('\ndemo_astar1 - apparent coordinates - jpl ephemeris\n');

fprintf('\n\nplease input a UTC calendar date\n');

[month, day, year] = getdate;

[uthr, utmin, utsec] = gettime;

day = day + uthr / 24 + utmin / 1440 + utsec / 86400;

jdutc = julian(month, day, year);

[filename, pathname] = uigetfile('*.dat', 'please select an input data file');

[fid, sname, ram, decm, pmra, pmdec, parlax, radvel] = readstar(filename);

% compute tdt julian date

jdtdt = utc2tdt(jdutc);

% define earth as central body

icent = 3;

while(1)
    
    fprintf('\ncoordinate type menu\n');
    
    fprintf('\n <1> geocentric\n');
    
    fprintf('\n <2> topocentric');
    
    fprintf('\n\nplease select coordinate type\n');
    
    ctype = input('? ');
    
    if (ctype == 1 || ctype == 2)
        break;
    end
    
end

if (ctype == 1)
    
    topo = 0;
    
else
    
    topo = 1;
    
end

if (topo == 1)
    
    [obslat, obslong, obsalt] = getobs;
    
    % define observer east longitude (degrees)
    
    glon = obslong * rtd;
    
    % define observer geodetic latitude (degrees)
    
    glat = obslat * rtd;
    
    % define observer altitude (meters)
    
    ht = 1000 * obsalt;
    
    [d, m, s, latstr] = deg2dms (glat);
    
    [d, m, s, longstr] = deg2dms (glon);
    
else
    
    glon = 0;
    
    glat = 0;
    
    ht = 0;
    
end

% calculate apparent coordinates

[ra, dec] = apstar1 (jdtdt, jdutc, icent, topo, glon, glat, ht, ...
    ram, decm, pmra, pmdec, parlax, radvel);

[h, m, s, rastr] = hrs2hms (ra);

[d, m, s, decstr] = deg2dms (dec);

% print results

if (topo == 0)
    
    fprintf ('\napparent geocentric coordinates \n\n');
    
else
    
    fprintf ('\napparent topocentric coordinates \n\n');
    
end

disp(sname);

[cdstr, utstr] = jd2str(jdutc);

fprintf('\nUTC calendar date        ');

disp(cdstr);

fprintf('\nUTC time                 ');

disp(utstr);

fprintf('\nUTC julian date        %14.4f \n', jdutc);

if (topo == 1)
    
    fprintf('\nobserver latitude        ');
    
    disp(latstr);
    
    fprintf('\nobserver east longitude  ');
    
    disp(longstr);
    
    fprintf('\nobserver altitude      %14.4f meters\n', 1000 * obsalt);
    
end

fprintf('\nright ascension          ');

disp(rastr);

fprintf('\ndeclination              ');

disp(decstr);

fprintf('\n');