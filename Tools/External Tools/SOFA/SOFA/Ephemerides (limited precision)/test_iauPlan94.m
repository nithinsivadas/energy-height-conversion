%  - - - - - - - - -
%   t _ p l a n 9 4
%  - - - - - - - - -
%
%  Test iauPlan94 function.
%
%  This revision:  2008 November 28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long g

pv = iauPlan94(2400000.5, -320000, 3)

pv = iauPlan94(2400000.5, 43999.9, 1)

