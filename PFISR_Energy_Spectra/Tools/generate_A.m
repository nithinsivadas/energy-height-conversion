function [A] = generate_A(h,E,latitude,longitude,time)
%generate_A
%   Calculates the A matrix: altitude profiles of  production rate per 
%   incident monoenergteic electron beams

% Functions required
% msis.m
% AIDA_TOOLS
% All of these are in Folder: PFISR_Energy_Spectra/Tools

% Input:
%  h    - altitude (km), double array [nh x 1]
%  E    - energy (eV) grid of the desired output spectra, double array
%         [nE x 1]
%  nO   - number density (m^-3) of atomic Oxygen [nh x 1]
%  nN2  - number density (m^-3) of molecular Nitrogen [nh x 1]
%  nO2  - number density (m^-3) of molecular Oxygen [nh x 1]
%  OPTS - Structure with options controlling the inversion
%         procedure. Current fields are:
%         OPTS.isotropic = 1; Run with isotropic electron precipitation,
%                     set to 2 for field-aligned and anything
%  Output:
% A     - is a matrix [nhxnE] which has the units (cm-1 eV)
%         Note that q = A*phi = [cm-1 eV]*[cm-2 s-1 eV-1]
% See also ISR2IeofE3 ionizationMatrixSemeter

%Default values to be replaced by switch case and input arguments

if nargin<5
    time        = datenum([2008 03 26 10 00 00]);
end
if nargin<4
    longitude   = -147.5; % Degrees.
end
if nargin<3
    latitude    = 65; % Degrees.
end




[D, T, F10_7_USED, AP_USED] = msis(time, latitude, longitude, h);
nO      = D(:,1);
nN2     = D(:,2);
nO2     = D(:,3);
density = D(:,4);

%% Estimating the Production Rates per Energy and Altitude using AIDA_TOOLS
[A] = ionization_profile_matrix(h,E,nO,nN2,nO2,density,1);


end

