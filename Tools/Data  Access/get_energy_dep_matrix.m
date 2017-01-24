function [A] = get_energy_dep_matrix(alt,energyBin,latitude,longitude,time)
% get_energy_dep_matrix.m
%   Calculates the energy deposition matrix matrix: altitude profiles of  production rate per 
%   incident monoenergteic electron beams

% Functions required
% msis.m
% AIDA_TOOLS: ionization_profile_matrix.m (and its following sub-functions)
% All of these are in Folder: PFISR_Energy_Spectra/Tools

% Input:
%  alt    		: altitude (km), double array [nh x 1]
%  energyBin    : energy (eV) grid of the desired output spectra, double array
%         		  [nE x 1]
%  latitude     : Latitude where the measurement is being made
%  longitude    : Longitude where the measurement is being made
%  time         : matlab units

%  Output:
%  A    		: is a matrix [nhxnE] which has the units (m-1 eV)
%         			Note that q = A*phi = [m-1 eV]*[m-2 s-1 eV-1]

%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------

	if nargin<5
	    time        = datenum([2008 03 26 10 00 00]);
	end
	if nargin<4
	    longitude   = -147.5; % Degrees. 
	end
	if nargin<3
	    latitude    = 65; % Degrees.
	end

	[D, T, F10_7_USED, AP_USED] = msis(time, latitude, longitude, alt);
	nO      = D(:,1);
	nN2     = D(:,2);
	nO2     = D(:,3);
	density = D(:,4);

	%% Estimating the Production Rates per Energy and Altitude using AIDA_TOOLS
	
	%  ionization_profile_matrix.m details
	%  nO   		: number density (m^-3) of atomic Oxygen [nh x 1]
	%  nN2  		: number density (m^-3) of molecular Nitrogen [nh x 1]
	%  nO2  		: number density (m^-3) of molecular Oxygen [nh x 1]
	%  OPTS 		: Structure with options controlling the inversion
	%         		procedure. Current fields are:
	%         		OPTS.isotropic = 1; Run with isotropic electron precipitation,
	%               set to 2 for field-aligned and anything

	[A] = ionization_profile_matrix(alt,energyBin,nO,nN2,nO2,density,1); % [cm-1 eV]
	
	% Converting to SI units
	A = (100)*A; %[m-1 eV]

end

