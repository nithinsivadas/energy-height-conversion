function [energyFlux] = num_to_energy(numberFlux,time,energyBins)

%num_to_energy_flux.m Converting differential number flux 
%to differential energy flux

% Input:
% numberFlux : NXM matrix [cm^-2 s^-1 eV^-1] or [m^-2 s^-1 eV^-1]
% time       : 1xM matrix [s]
% energyBins : Nx1 energy [eV]

% Output: 
% energyFlux : NxM matrix [eV cm^-2 sr^-1 s^-1 eV^-1] or [eV m^-2 sr^-1 s^-1 eV^-1]

%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

	[T,X]        = meshgrid(time,energyBins./(2*pi));
	energyFlux        = numberFlux.*X;

end

