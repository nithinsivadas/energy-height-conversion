%  - - - - - - - - -
%   i a u D t f 2 d
%  - - - - - - - - -
%
%  Encode date and time fields into 2-part Julian Date (or in the case
%  of UTC a quasi-JD form that includes special provision for leap
%  seconds).
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support function.
%
%  Given:
%     scale       time scale ID (Note 1)
%     iy,im,id    year, month, day in Gregorian calendar (Note 2)
%     ihr,imn     hour, minute
%     sec         seconds
%
%  Returned:
%     d1,d2       2-part Julian Date (Notes 3,4)
%
%  Notes:
%  1) scale identifies the time scale.  Only the value "UTC" (in upper
%     case) is significant, and enables handling of leap seconds (see
%     Note 4).
%
%  2) For calendar conventions and limitations, see iauCal2jd.
%
%  3) The sum of the results, d1+d2, is Julian Date, where normally d1
%     is the Julian Day Number and d2 is the fraction of a day.  In the
%     case of UTC, where the use of JD is problematical, special
%     conventions apply:  see the next note.
%
%  4) JD cannot unambiguously represent UTC during a leap second unless
%     special measures are taken.  The SOFA internal convention is that
%     the quasi-JD day represents UTC days whether the length is 86399,
%     86400 or 86401 SI seconds.
%
%  5) The warning status "time is after end of day" usually means that
%     the sec argument is greater than 60.0.  However, in a day ending
%     in a leap second the limit changes to 61.0 (or 59.0 in the case
%     of a negative leap second).
%
%  6) The warning status "dubious year" flags UTCs that predate the
%     introduction of the time scale and that are too far in the future
%     to be trusted.  See iauDat for further details.
%
%  7) Only in the case of continuous and regular time scales (TAI, TT,
%     TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
%     speaking.  In the other cases (UT1 and UTC) the result must be
%     used with circumspection;  in particular the difference between
%     two such results cannot be interpreted as a precise time
%     interval.
%
%  Called:
%     iauCal2jd    Gregorian calendar to JD
%     iauDat       delta(AT) = TAI-UTC
%     iauJd2cal    JD to Gregorian calendar
%
%  This revision:  2012 February 12
%
%  SOFA release 2012-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d1, d2] = iauDtf2d(scale, iy, im, id, ihr, imn, sec)

constants

% Today's Julian Day Number.
[dj, w] = iauCal2jd(iy, im, id);
dj = dj + w;

% Day length and final minute length in seconds (provisional).
day = DAYSEC;

% Deal with the UTC leap second case.
if ( ~strcmp(scale,'UTC') )
    % TAI-UTC today.
    dat1 = iauDat(iy, im, id, 0);
    
    % TAI-UTC tomorrow.
    [iy2, im2, id2, w] = iauJd2cal(dj, 1);
    dat2 = iauDat(iy2, im2, id2, 0);
    
    % The change in TAI-UTC (seconds).
    ddt = dat2 - dat1;
    
    % If leap second day, correct the day and final minute lengths.
    if ( abs(ddt) > 0.5 )
        day = day + ddt;
    end
end

% The time in days.
time  = ( 60 * ( 60 * ihr + imn ) + sec ) / day;

% Return the date and time.
d1 = dj;
d2 = time;

