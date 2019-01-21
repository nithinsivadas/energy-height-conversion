%  - - - - - - - -
%   t _ g d 2 g c
%  - - - - - - - -
%
%  Test iauGd2gc function.
%
%
%  Called:  iauGd2gc, viv, vvd
%
%  This revision:  2012 February 23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long g

e = 3.1;
p = -0.5;
h = 2500.0;
% WGS84
xyz = iauGd2gc(1, e, p, h)
% GRS80
xyz = iauGd2gc(2, e, p, h)
% WGS72
xyz = iauGd2gc(3, e, p, h)

