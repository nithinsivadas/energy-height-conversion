function [ numberFlux ] = kappa_j(x,energyBin)
%kappa_j.m Normalized flux (j(energyBin)) for particles with bulk energy
%(/velocity) meanEnergy,T,plasmaDensity, and kappa value k.

% Input:
% energyBin  	- Energy variable in [eV]
% meanEnergy 	- Mean energy in [eV]
% T  			- Temperature in [K]
% kappa  		- Kappa (unitless)
% plasmaDensity - density of plasma [m-3]
% m  			- mass of plasma [kg] Also 1 eV/(cm s-1)^2  = 1.6E-15 kg!  

% Output
% numberFlux  - Differential Number Flux [eV-1 cm-2 s-1]
%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------



meanEnergy=(x(1));
T=(x(2));
kappa=(x(3));
plasmaDensity=(x(4))*10^-6; %[cm^-3]
m=1;
m = m*(1.602*10^-15)^-1; %[eV/(cm s-1)^2]
numberFlux = (plasmaDensity*(energyBin).^0.5).*((2*pi*sqrt(2*m))^-1).*kappa_E(energyBin,meanEnergy,T,kappa);

% converting to SI units
numberFlux=numberFlux*10^4;

end

