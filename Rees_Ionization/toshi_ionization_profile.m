%% Generating altitude profiles of plasma density created by energetic precipitation from 3000 km altitude
% http://people.atmos.ucla.edu/ying/Spectrum/
% Isotropic, and within the polar cap in norther winter seasons

% Note that AIDA tool kit or the Sergienko and Ivanov 1993 model does not provide
% production rates for electron precipitation below 48.9 eV! (See the
% function ionization_profile_from_flux->energy_deg.m -> albedo.m or -> lambda.m


time        = datenum([2001 02 11 5 18 20]); % Must be a scalar or string.
latitude    = 66.95; % Degrees.
longitude   = -50.94; % Degrees.
h           = 80:5:500; % km. Must be a scalar.
%M           = importdata('Semeter2005_100eV.dat',' ',1);
%E1          = M.data(end:-1:1,1);
%dE          = diff(E1);
%E           = E1(2:end); %[eV cm-2 s-1 sr-1 eV-1]
%flux        = 2*pi*dE.*M.data(end:-1:2,2)./E; % [cm-2 s-1] Converting energy flux into number flux
% E           = [98,99,100,101,102]';
E           = (998:1:1002)';
flux        = [0,0,1.0e8,0,0]';

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

[A] = ionization_profile_matrix(h',E,nO,nN2,nO2,density,1);

% OPTS=1 : Isotropic precipitation

%% Estimating the Ionization profile
q=nansum(A(:,:).*repmat(flux,size(h))',2);

alpha_fit = (2.5*10^-12)*exp(-h/51.2); % [m^3/s] Effective Recombination Coefficient
alpha_fit = alpha_fit*10^6;           % [cm^3/s]

% Assuming q = alpha_fit*n_e^2
n_e = (q./alpha_fit').^0.5;            % [cm^-3]

%% Plotting the Electron Density Profile
figure;
semilogx(q*10^6,h);
xlabel('Production Rate [m^-^3s^-1]')
ylabel('Altitude [km]');
title ('2001-02-11 05:18:20 UT');
grid on;
xlim([10^7 10^11]);
ylim([50 300]);

figure;
semilogx(density,h);
xlabel('Density [m^-^3]')
ylabel('Altitude [km]');
title ('2001-02-11 05:18:20 UT');
grid on;
xlim([10^-12 10^-4]);
ylim([50 300]);


%% Writing the data on a text file
str='production_rate_profile_Semeter_2005.txt';
str1='#Altitude [km] Electron Density [cm-3]';
dlmwrite(str,str1,'-append', 'delimiter', '','newline', 'pc'); 
str2=[(h'),(n_e)];
dlmwrite(str,str2,'-append', 'delimiter', '\t','newline', 'pc'); 

