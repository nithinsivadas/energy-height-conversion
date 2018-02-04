function [range,az,el] = cart2sphere(x,y,z,varargin)
% sphere2cart.m
% by John Swoboda
% This is a simple transform that will change cartisian coordinates, in
% to spherical in radians or in degrees depending on the input type. It is assumed that the elevation is measured from
% the z=0 plane and not the x=y=0 line.
% Inputs 
% x,y,z - Arrays that hold the cartisian coordinates. Must all be the same
% size.
% type - Optional string argument that is either 'deg' or 'rad'. The
% default is out put in radians
% Outputs
% range - An array of range points in whatever unit x, y and z will be in.
% az,el - Arrays that hold the az and el components both in radians or degrees 
% depending on the input type. Elevation will be referenced from z=0 plane. 
% Will be same size as range.
% Outputs
if nargin==3
    r2d = 1;
elseif strcmpi(varargin{1},'deg')
    r2d = 180/pi;
elseif strcmpi(varargin{1},'rad')
    r2d=1;
end

range = sqrt(x.^2+y.^2+z.^2);
az = atan2(x,y)*r2d;
el = asin(z./range)*r2d;