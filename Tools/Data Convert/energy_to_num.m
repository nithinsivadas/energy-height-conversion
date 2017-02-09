function [numberFlux] = energy_to_num(energyFlux, time, energyBin)
%% energy_to_num.m Converting differential energy flux to differential number flux

% Input:
% energyFlux : NxM matrix [eV cm^-2 sr^-1 s^-1 eV^-1] or [eV m^-2 sr^-1 s^-1 eV^-1]
% time       : 1xM matrix [s]
% energyBins : Nx1 energy [eV]

% Output: 
% numberFlux : NXM matrix [cm^-2 s^-1 eV^-1] or [m^-2 s^-1 eV^-1]
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

	[T,X]        = meshgrid(time,(energyBin./(2*pi)).^-1);
	numberFlux   = energyFlux.*X;

	    [isThereNAN, totalNAN] = check_nan(numberFlux);
end

