function [fid, sname, ram, decm, pmra, pmdec, parlax, radvel] = readstar(filename)

% read star data file

% required by star*.m

% input

%  filename = name of star data file

% output

%  sname  = star name
%  ram    = J2000 right ascension (hours)
%  decm   = J2000 declination (degrees)
%  pmra   = J2000 proper motion in right ascension
%           (seconds/Julian century)
%  pmdec  = J2000 proper motion in declination
%           (arcseconds/Julian century)
%  parlax = parallax (arcseconds)
%  radvel = radial velocity (kilometers/second)
%  fid    = file id

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    clc; home;
    fprintf('\n\n  error: cannot find this file!!');
    keycheck;
    return;
end

% read 20 lines of data file

for i = 1:1:20
    cline = fgetl(fid);
    switch i
        case 2
            sname = cline;
        case 5
            ram = str2double(cline);
        case 8
            decm = str2double(cline);
        case 11
            pmra = str2double(cline);
        case 14
            pmdec = str2double(cline);
        case 17
            parlax = str2double(cline);
        case 20
            radvel = str2double(cline);
    end
end

status = fclose(fid);

