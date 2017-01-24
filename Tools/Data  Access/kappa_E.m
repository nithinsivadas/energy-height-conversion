function [ noParticles ] = kappa_E(energyBin, meanEnergy, T, kappa)
% kappa_E.m Normalized Energy distribution function (noParticles) for particles with mean energy
% E0,T, and kappa value kappa.

% Input:
% energyBin  - Energy variable in [eV]
% meanEnergy - Mean energy in [eV]
% T  - Temperature in [kappa]
% kappa  - Kappa (unitless)

% Output:
% noParticles - no. of particles per Energy [eV-1]
%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------

eV = 1.60218E-19;
kb = 1.38*10^-23; %Boltzmann Constant
energyBin  = energyBin*eV;
meanEnergy = meanEnergy*eV;

a1 = pi*meanEnergy*kb*T*(kappa-1.5);
a2 = gamma(kappa)/gamma(kappa-0.5);
a3 = 1+((energyBin.^0.5 - meanEnergy^0.5).^2)/((kappa-1.5)*kb*T);
a4 = 1+((energyBin.^0.5 + meanEnergy^0.5).^2)/((kappa-1.5)*kb*T);

% No. of particles per Joule of Energy
noParticles = 0.5*(a1^-0.5)*a2*(a3.^-kappa - a4.^-kappa);  % [J^-1]

% No. of particles per eV
noParticles = noParticles*eV;                % [eV^-1]

end

