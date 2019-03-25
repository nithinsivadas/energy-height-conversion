%% Relativisitc maxwellian fitting
E = 30:1:1000; % keV
T = 150; %keV
n = 4e-4; %cm-3

f1 = rel_maxwellian(E,5e-3,30);
f2 = rel_maxwellian(E,1e-4,170);

figure; 

loglog(E,f1,'r');
hold on;

loglog(E,f2,'b');
hold on;
loglog(E,f1+f2,'k');
legend('Ts = 25 keV', 'Th = 170 keV', 'Total');
set(gca,'XTick',[30, 100, 200, 300, 1000]);
xlabel('Energy (E) in KeV');
ylabel('f(E) - Relativistic maxwellian distribution');
grid on;

figure; 
loglog(E,100*f1./(f1+f2),'r');
hold on;
loglog(E,100*f2./(f1+f2),'b');
title('% of total flux');
legend('Ts = 25 keV', 'Th = 170 keV');
set(gca,'XTick',[30, 100, 200, 300, 1000],'YTick',[1,10,30,50,100]);
xlabel('Energy (E) in KeV');
ylabel('Percentage');
ylim([1 100]);
grid on;

function f = rel_maxwellian(E, n, T)
C = define_universal_constants();

% Convert to SI units
E = E.*1e+3.*C.e; % J 
T = T.*1e+3.*C.e./C.kb; % keV to K;
n = n*1e+6;
m = C.me;
a = (m.*C.c^2)./(C.kb.*T);
k2 = besselk(2,a);
f = n.*((2*m.*C.c.^2).^0.5)./(4*pi.*(C.kb.*T).^2);
f = f*(exp(-E./(C.kb.*T)))./(a.*exp(a).*k2);


end