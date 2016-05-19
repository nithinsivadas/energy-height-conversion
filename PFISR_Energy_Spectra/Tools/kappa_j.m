function [ j ] = kappa_j(Eb,T,k,n,m,E)
%kappa_j.m Normalized flux (j(E)) for particles with bulk energy
%(/velocity) Eb,T,n, and kappa value k.

% Input:
% E  - Energy variable in [eV]
% Eb - Mean energy in [eV]
% T  - Temperature in [K]
% k  - Kappa (unitless)
% n  - density of plasma [cm-3]
% m  - mass of plasma [kg] Also 1 eV/(cm s-1)^2  = 1.6E-15 kg!  

% Output
% j  - Differential Number Flux [eV-1 cm-2 s-1]

m = m*(1.602*10^-15)^-1; %[eV/(cm s-1)^2]
j = (n*(E).^0.5).*((2*pi*sqrt(2*m))^-1).*kappa_E(E,Eb,T,k);

end

