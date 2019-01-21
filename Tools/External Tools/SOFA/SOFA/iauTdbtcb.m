%  - - - - - - - - - -
%   i a u T d b t c b
%  - - - - - - - - - -
%
%  Time scale transformation:  Barycentric Dynamical Time, TDB, to
%  Barycentric Coordinate Time, TCB.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  canonical.
%
%  Given:
%     tdb1,tdb2      TDB as a 2-part Julian Date
%
%  Returned:
%     tcb1,tcb2      TCB as a 2-part Julian Date
%
%  Notes:
%  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
%     between the two arguments, for example where tdb1 is the Julian
%     Day Number and tdb2 is the fraction of a day.  The returned
%     tcb1,tcb2 follow suit.
%
%  2) The 2006 IAU General Assembly introduced a conventional linear
%     transformation between TDB and TCB.  This transformation
%     compensates for the drift between TCB and terrestrial time TT,
%     and keeps TDB approximately centered on TT.  Because the
%     relationship between TT and TCB depends on the adopted solar
%     system ephemeris, the degree of alignment between TDB and TT over
%     long intervals will vary according to which ephemeris is used.
%     Former definitions of TDB attempted to avoid this problem by
%     stipulating that TDB and TT should differ only by periodic
%     effects.  This is a good description of the nature of the
%     relationship but eluded precise mathematical formulation.  The
%     conventional linear relationship adopted in 2006 sidestepped
%     these difficulties whilst delivering a TDB that in practice was
%     consistent with values before that date.
%
%  3) TDB is essentially the same as Teph, the time argument for the
%     JPL solar system ephemerides.
%
%  Reference:
%     IAU 2006 Resolution B3
%
%  This revision:  2011 May 14
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tcb1, tcb2] = iauTdbtcb(tdb1, tdb2)

constants

% 1977 Jan 1 00:00:32.184 TT, as two-part JD
t77td = DJM0 + DJM77;
t77tf = TTMTAI/DAYSEC;

% TDB (days) at TAI 1977 Jan 1.0
tdb0 = TDB0/DAYSEC;

% TDB to TCB rate
elbb = ELB/(1-ELB);

% Result, preserving date format but safeguarding precision.
if (tdb1 > tdb2)
    d = t77td - tdb1;
    f  = tdb2 - tdb0;
    tcb1 = tdb1;
    tcb2 = f - ( d - ( f - t77tf ) ) * elbb;
else
    d = t77td - tdb2;
    f  = tdb1 - tdb0;
    tcb1 = f + ( d - ( f - t77tf ) ) * elbb;
    tcb2 = tdb2;
end

