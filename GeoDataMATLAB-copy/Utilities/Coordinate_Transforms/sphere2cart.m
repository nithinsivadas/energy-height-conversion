function [x,y,z] = sphere2cart(range,az,el)
% sphere2cart.m
% by John Swoboda
% This is a simple transform that will change sphereical coordinates, in
% radians to cartisian. It is assumed that the elevation is measured from
% the z=0 plane and not the x=y=0 line.
% Inputs 
% range - An array of range points in whatever unit x, y and z will be in.
% az,el - Arrays that hold the az and el components both in radians and
% elevation will be referenced from z=0 plane. Must be same size as range
% Outputs
% x,y,z - Arrays the same size of range, az and el. Will be whatever units
% range was in.
kx = sin(az) .* cos(el);
ky = cos(az) .* cos(el);
kz = sin(el);

x = range.*kx;
y = range.*ky;
z = range.*kz;