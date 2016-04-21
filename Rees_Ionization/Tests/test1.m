%% Test1.m
% This routine produces the production rate altitude profile for an 
% incident electron beam of 49 keV of flux 1.0e8 cm-2 s-1. It is to 
% be verified with figure 2 of Semeter 2005. 

%% What we learnt from this test
% The flux value to be multiplied with [A] to get the production rate q
% ought to be divided by the energy bin. i.e flux/dE [cm-2 s-1 eV-1]

%% Initialization

clear all;

addpath('/home/nithin/Documents/git-repos/energy-height-conversion/Rees_Ionization')

time        = datenum([2001 02 11 5 18 20]); 	% 
latitude    = 67; 								% [Deg]
longitude   = -50; 							    % [Deg]
h           = 80:0.5:500; 						% [km] 
E           = (48998:1:49002);                    
dE          = diff(E);
dE(1)       = dE(1);
dE(end+1)   = dE(end);
flux        = [0,0,1.0e8,0,0]./dE;              % [cm-2 s-1 eV-1]

%% MSIS-90 Model
[D, T, F10_7_USED, AP_USED] = msis(time, latitude, longitude, h);

nO      = D(:,1);
nN2     = D(:,2);
nO2     = D(:,3);
density = D(:,4);

%       1. O: Atomic oxygen (O) number density in m^-3.
%       2. N2: Molecular nitrogen (N_2) number density in m^-3.
%       3. O2: Molecular oxygen (O_2) number density in m^-3.
%       4. MASS_DENSITY: Total mas density in kg/m^3.

%% Estimating the Production Rates per Energy and Altitude using AIDA_TOOLS

[A] = ionization_profile_matrix(h',E',nO,nN2,nO2,density,1);
% OPTS=1 : Isotropic precipitation

%% Estimating the Ionization profile
q=A*flux'; %cm^-3 s^-1

%% Plotting the Production Rates
figure;
semilogx(q*10^6,h);
xlabel('Production Rates [m^-^3 s^-^1]')
ylabel('Altitude [km]');
title ('2005-03-11 05:18:20 UT');
grid on;
xlim([10^7 10^11]);
