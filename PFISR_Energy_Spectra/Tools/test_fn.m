function [ y ] = test_fn( x,xdata )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
m = (9.11*10^-31);
E = xdata;
Eb=x(1);
T =x(2);
k =x(3);
n =3*10^3;
y=kappa_j(Eb,T,k,n,m,E);

end

