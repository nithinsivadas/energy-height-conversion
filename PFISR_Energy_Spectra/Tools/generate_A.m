function [A] = generate_A(h,E,latitude,longitude,time)
%generate_A
%   Calculates the A matrix: altitude profiles of  production rate per 
%   incident monoenergteic electron beams

% Functions required
% msis.m
% AIDA_TOOLS
% All of these are in Folder: PFISR_Energy_Spectra/Tools

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

