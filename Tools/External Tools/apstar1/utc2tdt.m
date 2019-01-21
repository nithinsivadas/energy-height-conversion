function jdtdt = utc2tdt (jdutc)

% convert UTC julian date to TDT julian date

% input

%  jdutc   = UTC julian date

% output

%  jdtdt = TDT julian date 

% Reference Frames in Astronomy and Geophysics
% J. Kovalevsky et al., 1989, pp. 439-442

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find current number of leap seconds

leapsecond = findleap(jdutc);

% TDT julian date
    
corr = (leapsecond + 32.184) / 86400.0;
    
jdtdt = jdutc + corr;
    
