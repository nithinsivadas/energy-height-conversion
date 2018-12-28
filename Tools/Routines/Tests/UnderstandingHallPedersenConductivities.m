%% Understanding the Hall and Pedersen Conductivities

alt = logspace(log10(50),log10(600),100);
Ne  = ones(size(alt)) * 10^10;

[data, input] = get_conductivity_v2( alt, Ne,...
    [], [], [], 0, {'all'}, true);

%%
k_i = interp1(data.altitude',data.ionMobility,alt');
k_e = interp1(data.altitude',data.electronMobility,alt');
C = interp1(data.altitude,input.ionConcentration,alt);
sigmaPi = zeros(length(alt),1);
for i=1:1:size(k_i,2)
    sigmaPi = sigmaPi + C(:,i).*k_i(:,i)./(1+(k_i(:,i).^2));
end
sigmaPe = (k_e)./(1+k_e.^2);

sigmaHi = zeros(length(alt),1);
for i=1:1:size(k_i,2)
    sigmaHi = sigmaHi + C(:,i).*(k_i(:,i).^2)./(1+(k_i(:,i).^2));
end
sigmaHe = (k_e.^2)./(1+k_e.^2);

%%
figure;
semilogx(sigmaPe,alt,'.-b'); hold on;
semilogx(sigmaPi,alt,'*-r'); hold on;
semilogx(real(sigmaPi+sigmaPe),alt,'--k');
legend('e^-','i^+','Total');
xlim([1e-3 1]);
ylim([60 600]);
title({'Normalized contribution from electrons and ions','to \sigma_P'});
set(gca,'XTick',[1e-3,1e-2,1e-1,1],'XTickLabel',{'0.001','0.01','0.1','1'});
ylabel('[km]');
grid on;

figure;
semilogx(sigmaHe,alt,'.-b'); hold on;
semilogx(sigmaHi,alt,'*-r'); hold on;
semilogx(real(-sigmaHi+sigmaHe),alt,'--k');
legend('e^-','i^+','Total');
xlim([1e-3 2]);
ylim([60 600]);
title({'Normalized contribution from electrons and ions','to \sigma_H'});
set(gca,'XTick',[1e-3,1e-2,1e-1,1],'XTickLabel',{'0.001','0.01','0.1','1'});
ylabel('[km]');
grid on;