clear all;

energyBin = logspace(3,6,100)';
n = length(energyBin);
alt  = (60:1:150)';
lat  = 65;
long = -147.5;
time = datenum([2008 03 26 10 00 00]);

A = get_energy_dep_matrix(alt, energyBin, repmat(lat,length(alt),1), repmat(long,length(alt),1), time);

numFlux = diag(zeros(1,length(energyBin)));
numFlux(find_altitude(energyBin,100000)) = 10^8;

q = A*numFlux;

%% GLOW
time = datenum('2008-03-26 11:00');
[f107a, f107, aph] = f107_aph(time);
Ap = aph(:,1);
f107p = (f107a + f107)./2;
iono = glowenergy(time, lat, long, f107a, f107, f107p, Ap, energyBin, numFlux(:,1)*1e-4);

%%
figure; 
plot(energyBin,numFlux(:,1)+1);
set(gca,'XScale', 'log','YScale','log');
xlabel('Energy [eV]');
ylabel ('Electron number flux [m^{-2} s^{-1} eV^{-1}]');

figure; 
plot(q(:,1), alt); xlim([10^6, 10^12]); 
set(gca,'XScale','log');
xlabel('Prodcution rate [cm^{-2}s^{-1}]');
ylabel('Altitude [km]');

figure;
plot(iono.A4278,iono.alt,'b');
hold on;
plot(iono.A5577,iono.alt,'g');
hold on;
plot(iono.A6300,iono.alt,'r');
set(gca,'XScale','log');
ylim([alt(1),alt(end)]);
xlabel('Volume emission rates cm^{-3} s^{-1}');
ylabel('Altitude [km]');
legend('A4278','A5577','A6300');

figure;
plot(iono.pedersen,iono.alt,'r');
hold on;
plot(iono.hall,iono.alt,'k');
set(gca,'XScale','log');
ylim([alt(1),alt(end)]);
xlabel('Conductivity [S]');
ylabel('Altitude [km]');
legend('\sigma_P','\sigma_H');