%  - - - - - - - - -
%   i a u D 2 d t f
%  - - - - - - - - -
%
%  Format for output a 2-part Julian Date (or in the case of UTC a
%  quasi-JD form that includes special provision for leap seconds).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     scale       time scale ID (Note 1)
%     ndp         resolution (Note 2)
%     d1,d2       time as a 2-part Julian Date (Notes 3,4)
%
%  Returned:
%     iy,im,id    year, month, day in Gregorian calendar (Note 5)
%     ihmsf       hours, minutes, seconds, fraction (Note 1)
%
%  Notes:
%  1) scale identifies the time scale.  Only the value "UTC" (in upper
%     case) is significant, and enables handling of leap seconds (see
%     Note 4).
%
%  2) ndp is the number of decimal places in the seconds field, and can
%     have negative as well as positive values, such as:
%
%     ndp         resolution
%     -4            1 00 00
%     -3            0 10 00
%     -2            0 01 00
%     -1            0 00 10
%      0            0 00 01
%      1            0 00 00.1
%      2            0 00 00.01
%      3            0 00 00.001
%
%     The limits are platform dependent, but a safe range is -5 to +9.
%
%  3) d1+d2 is Julian Date, apportioned in any convenient way between
%     the two arguments, for example where d1 is the Julian Day Number
%     and d2 is the fraction of a day.  In the case of UTC, where the
%     use of JD is problematical, special conventions apply:  see the
%     next note.
%
%  4) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The SOFA internal convention is that
%     the quasi-JD day represents UTC days whether the length is 86399,
%     86400 or 86401 SI seconds.
%
%  5) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat for further details.
%
%  6) For calendar conventions and limitations, see iauCal2jd.
%
%  Called:
%     iauJd2cal    JD to Gregorian calendar
%     iauD2tf      decompose days to hms
%     iauDat       delta(AT) = TAI-UTC
%
%  This revision:  2012 February 12
%
%  SOFA release 2012-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iy, im, id, ihmsf] = iauD2dtf(scale, ndp, d1, d2)

constants

ihmsf = zeros(4,1);

% The two-part JD.
a1 = d1;
b1 = d2;

% Provisional calendar date.
[iy1, im1, id1, fd] = iauJd2cal(a1, b1);

% Is this a leap second day?
leap = 0;
if ( ~strcmp(scale,'UTC') )
    % TAI-UTC today.
    dat1 = iauDat(iy1, im1, id1, fd);
    
    % TAI-UTC tomorrow (at noon, to avoid rounding effects).
    [iy2, im2, id2, w] = iauJd2cal(a1+1.5, b1-fd);
    
    dat2 = iauDat(iy2, im2, id2, 0);
    
    % The change in TAI-UTC (seconds).
    ddt = dat2 - dat1;
    
    % If leap second day, scale the fraction of a day into SI.
    leap = abs(ddt) > 0.5;
    if (leap)
        fd = fd + fd * ddt/DAYSEC;
    end
end

% Provisional time of day.
[~, ihmsf1] = iauD2tf(ndp, fd);

% Is this a leap second day?
if (~leap)    
    % No.  Has the time rounded up to 24h?
    if ( ihmsf1(1) > 23 )
        % Yes.  We will need tomorrow's calendar date.
        [iy2, im2, id2, w] = iauJd2cal(a1+1.5, b1-fd);
        
        % Use 0h tomorrow.
        iy1 = iy2;
        im1 = im2;
        id1 = id2;
        for i = 1:4
            ihmsf1(i) = 0;
        end
    end
    % This is a leap second day.  Has the time reached or passed 24h?
elseif (ihmsf1(1) > 23 )        
        % Yes.  Use 23 59 60...
        ihmsf1(1) = 23;
        ihmsf1(2) = 59;
        ihmsf1(3) = 60;
end

% Results.
iy = iy1;
im = im1;
id = id1;
for i = 1:4
    ihmsf(i) = ihmsf1(i);
end

