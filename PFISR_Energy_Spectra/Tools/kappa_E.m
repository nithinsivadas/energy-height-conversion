function [ P_E ] = kappa_E(E, Eb, T, k)
%kappa_E.m Normalized Energy distribution function (P_E) for particles with mean energy
%E0,T, and kappa value k.

% Input:
% E  - Energy variable in [eV]
% E0 - Mean energy in [eV]
% T  - Temperature in [K]
% k  - Kappa (unitless)

% Output:
% P_E - no. of particles per Energy [eV-1]

eV = 1.60218E-19;
kb = 1.38*10^-23; %Boltzmann Constant
E  = E*eV;
Eb = Eb*eV;

a1 = pi*Eb*kb*T*(k-1.5);
a2 = gamma(k)/gamma(k-0.5);
a3 = 1+((E.^0.5 - Eb^0.5).^2)/((k-1.5)*kb*T);
a4 = 1+((E.^0.5 + Eb^0.5).^2)/((k-1.5)*kb*T);

% No. of particles per Joule of Energy
P_E = 0.5*(a1^-0.5)*a2*(a3.^-k - a4.^-k);  % [J^-1]

% No. of particles per eV
P_E = P_E*eV;                % [eV^-1]

end

