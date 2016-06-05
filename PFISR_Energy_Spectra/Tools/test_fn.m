function [ y ] = test_fn( x,xdata )
%test_fn A function that generates a kappa function
%   Detailed explanation goes here


% Input:
% E  - Energy variable in [eV]
% Eb - Mean energy in [eV]
% T  - Temperature in [K]
% k  - Kappa (unitless)
% n  - density of plasma [cm-3]
% m  - mass of plasma [kg] Also 1 eV/(cm s-1)^2  = 1.6E-15 kg!  

m = (9.11*10^-31);
E = xdata;
Eb=x(1);
T =x(2);
k =x(3);
n =3*10^3;
y=kappa_j([Eb,T,k,n],E);

end

