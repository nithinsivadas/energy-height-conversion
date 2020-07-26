%% Understanding the Hall and Pedersen Conductivities

alt = logspace(log10(50),log10(600),100);
Ne  = ones(size(alt)) * 10^10;

[data, input] = get_conductivity_v2( alt, Ne,...
    [], [], [], 0, {'all'}, true);

%%
NnNi = 0.5; % ratio of negative ions to positive ions
dmobility = 0.1; % The reduction in mobility due to heavier negative ions

k_i = interp1(data.altitude',data.ionMobility,alt');
k_e = interp1(data.altitude',data.electronMobility,alt');
C = interp1(data.altitude,input.ionConcentration,alt);
Cn = zeros(size(alt))';
Cn(alt<80) = NnNi;
nNe = 1-Cn; % concentration of electrons

sigmaPi = zeros(length(alt),1);
for i=1:1:size(k_i,2)
    sigmaPi = sigmaPi + C(:,i).*k_i(:,i)./(1+(k_i(:,i).^2));
end

sigmaPn = Cn.*(k_i(:,2).*dmobility)./(1+((k_i(:,2).^2).*dmobility));
sigmaPe = (k_e)./(1+k_e.^2);

sigmaHi = zeros(length(alt),1);
for i=1:1:size(k_i,2)
    sigmaHi = sigmaHi + C(:,i).*(k_i(:,i).^2)./(1+(k_i(:,i).^2));
end
sigmaHn = Cn.*((k_i(:,2).*dmobility).^2)./(1+((k_i(:,2).*dmobility).^2));
sigmaHe = (k_e.^2)./(1+k_e.^2);

%
figure;
semilogx(sigmaPe,alt,'.-b'); hold on;
semilogx(sigmaPi,alt,'*-r'); hold on;
semilogx(sigmaPn,alt,'--g'); hold on;
semilogx(real(sigmaPi./nNe+sigmaPe+Cn.*sigmaPn./nNe),alt,'-k');
% semilogx(real(sigmaPi+sigmaPe),alt,'--k');
legend('e^-','i^+','i^-','Total','Location','northwest');
xlim([1e-6 1]);
ylim([60 150]);
title({'Normalized contribution from electrons and ions','to \sigma_P'});
set(gca,'XTick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1],'XTickLabel',{'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'});
ylabel('[km]');
grid on;

figure;
semilogx(sigmaHe,alt,'.-b'); hold on;
semilogx(sigmaHi,alt,'*-r'); hold on;
semilogx(sigmaHn,alt,'--g'); hold on;
semilogx(real(-sigmaHi./nNe+sigmaHe+Cn.*sigmaHn./nNe),alt,'-k');
% semilogx(real(-sigmaHi+sigmaHe),alt,'--k');
legend('e^-','i^+','i^-','Total','Location','northwest');
xlim([1e-6 2]);
ylim([60 150]);
title({'Normalized contribution from electrons and ions','to \sigma_H'});
set(gca,'XTick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1],'XTickLabel',{'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'});
ylabel('[km]');
grid on;