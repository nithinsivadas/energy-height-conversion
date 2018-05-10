function [A] = get_energy_dep_matrix(alt,energyBin,latitude,longitude,time)
% get_energy_dep_matrix.m Calculates the energy deposition matrix : 
% altitude profiles of  production rate per incident monoenergteic electron beams

% Functions required
% 1. msis.m
% 2. AIDA_TOOLS: ionization_profile_matrix.m (and its following sub-functions)
% All of these are in Folder: PFISR_Energy_Spectra/Tools

%-------------------------------------------------------------------------
% Input:
%-------
%  alt    		: altitude (km), double array [nh x 1]
%  energyBin    : energy (eV) grid of the desired output spectra, double array
%         		  [nE x 1]
%  latitude     : Latitude where the measurement is being made
%  longitude    : Longitude where the measurement is being made
%  time         : matlab units
%-------------------------------------------------------------------------
%  Output:
%---------
%  A    		: is a matrix [nhxnE] which has the units (m-1 eV)
%         	      Note that q = A*phi = [m-1 eV]*[m-2 s-1 eV-1]
%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
    nAlt = length(alt);
	if nargin<5
	    time        = datenum([2008 03 26 10 00 00]);
	end
	if nargin<4
	    longitude   = -147.5*ones(length(alt),1); % Degrees. 
	end
	if nargin<3
	    latitude    = 65*ones(length(alt),1); % Degrees.
    end
    
    [F107A, F107, APH] = f107_aph(time);
    msis00Data = onera_desp_lib_msis('nrlmsise00',repmat(time,nAlt,1),...
        [alt,latitude,longitude],0,repmat(F107A,nAlt,1),...
        repmat(F107,nAlt,1),repmat(APH,nAlt,1));
    nO = msis00Data.O*10^6; %m^-3
    nN2 = msis00Data.N2*10^6; %m^-3
    nO2 = msis00Data.O2*10^6; %m^-3
    density = msis00Data.TotalMass*10^3; % kg m^-3
    
% 	[D, T, F10_7_USED, AP_USED] = msis(time, latitude, longitude, alt);
% 	nO      = D(:,1);
% 	nN2     = D(:,2);
% 	nO2     = D(:,3);
% 	density = D(:,4);

	%% Estimating the Production Rates per Energy and Altitude using AIDA_TOOLS
	
	%  ionization_profile_matrix.m details
	%  nO   		: number density (m^-3) of atomic Oxygen [nh x 1]
	%  nN2  		: number density (m^-3) of molecular Nitrogen [nh x 1]
	%  nO2  		: number density (m^-3) of molecular Oxygen [nh x 1]
    %  density      : mass density kg/m^3
	%  OPTS 		: Structure with options controlling the inversion
	%         		procedure. Current fields are:
	%         		OPTS.isotropic = 1; Run with isotropic electron precipitation,
	%               set to 2 for field-aligned and anything

	[A] = ionization_profile_matrix(alt,energyBin,nO,nN2,nO2,density,1); % [cm-1 eV]
	
	% Converting to SI units
	A = (100)*A; %[m-1 eV]

end

