%% Plot the peak energy per altitude
clear all;

energyBin = logspace(3,6,100)';
alt  = (60:1:150)';
lat  = 65;
long = -147.5;
time = datenum([2008 03 26 10 00 00]);

A = get_energy_dep_matrix(alt, energyBin, lat, long, time);

numFlux = diag(ones(1,length(energyBin)))*10^12;

q = A*numFlux;

[q_max, ialt] = max(q);

figure;
semilogx((q),alt); xlim([10^10 10^18]); hold on;
figure;
semilogx(energyBin,alt(ialt));