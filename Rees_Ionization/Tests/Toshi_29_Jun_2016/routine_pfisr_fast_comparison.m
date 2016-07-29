%% Comparing FAST & PFISR measurements of soft precipitation

time      = datenum([2008 03 04 06 35 20]);
latitude  = 66.95;  % Degrees
longitude =-50.93;  % Degrees

load pfisrfastcomparison/20080304_0636.mat 
load pfisrfastcomparison/20080304_0636_alpha_toshi.mat

%% Variales loaded: 
% h    - altitude          [km]
% ne   - electron density  [m-3]
% E    - energy bins       [eV]
% flux - energy flux       [eV/ cm2 s sr eV]
% alpha_toshi - effective recombination rate [cm3/s]
% h_alpha     - altitudes where alpha_toshi is calculated [km]
%% Calculating differential number flux
E    = fliplr(E(1:1:end));
flux = fliplr(flux(1:1:end));

dE = diff(E);
dE(1) = dE(1);
dE(end+1)   = dE(end);
flux        = flux./dE; % [cm-2 sr-1 eV-1]
flux        = flux*2*pi; % [cm-2 eV-1]

%% MSIS-90 Model
[D, T, F10_7_USED, AP_USED] = msis(time, latitude, longitude, h);

nO      = D(:,1);
nN2     = D(:,2);
nO2     = D(:,3);
density = D(:,4);

%% Estimating the Production Rates per Energy and Altitude using AIDA_TOOLS

[A] = ionization_profile_matrix(h,E,nO,nN2,nO2,density,1);

%% Estimating the Ionization profile
q_est=nansum(A(:,:).*repmat(flux',size(h)),2); % [cm-3 s-1]
q_est=q_est*10^6;

% alpha_fit = (2.5*10^-12)*exp(-h/51.2); % [m^3/s] Effective Recombination Coefficient
alpha_fit = interp1(h_alpha,alpha_toshi*10^-6,h);

% alpha_fit = alpha_fit*10^6;           % [cm^3/s]

% Assuming q = alpha_fit*n_e^2
ne_est = (q_est./alpha_fit).^0.5;      % [m^-3]


%% Plotting the estimate ne and measured ne
figure;
semilogx(ne_est,h,'-r');
hold on;
semilogx(ne,h,'-black');
legend('Estimated ne','Measured ne');
xlabel('Density [m^-^3]')
ylabel('Altitude [km]');
title ('2008-03-04 06:36 UT');
grid on;
% xlim([10^5 10^11]);
ylim([100 700]);
