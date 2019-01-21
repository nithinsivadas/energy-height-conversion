%  - - - - - - - -
%   t _ e p v 0 0
%  - - - - - - - -
%
%  Test iauEpv00 function.
%
%  This revision:  2008 November 28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long g

[pvh, pvb] = iauEpv00(2400000.5, 53411.52501161)

