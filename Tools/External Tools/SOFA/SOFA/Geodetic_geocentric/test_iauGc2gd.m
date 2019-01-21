%  - - - - - - - -
%   t _ g c 2 g d
%  - - - - - - - -
%
%  Test iauGc2gd function.
%
%  Called:  iauGc2gd
%
%  This revision:  2012 February 23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long g

xyz = [2e6, 3e6, 5.244e6];
% WGS84
[e, p, h] = iauGc2gd(1, xyz)
% GRS80
[e, p, h] = iauGc2gd(2, xyz)
% WGS72
[e, p, h] = iauGc2gd(3, xyz)

