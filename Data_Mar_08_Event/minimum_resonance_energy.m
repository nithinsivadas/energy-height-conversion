%% Minumum resonant energy calculation
% Script to calculate and plot the minimum energy of electrons that may
% resonate with EMIC waves

% Uncertain parameters 
% amp : the wave power for maximum cut-off
% fraction: composition of the ions in the plasmasheet
clear all;

time=linspace(datenum('2008-03-26/08:00'),datenum('2008-03-26/13:00'),1000);
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Bfield_thd.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/thd_density.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Waves/GSM_xwv.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Waves/FGM_xwv.mat');


B   = interp1(thd_FGM.time,thd_FGM.Btot,time,'linear','extrap'); % nT
ne  = interp1(thd_density.time,thd_density.e,time,'linear','extrap'); % cm-3
ni  = interp1(thd_density.time,thd_density.i,time,'linear','extrap'); % cm-3
wvpow = interp1(xwvfgm.time,xwvfgm.pow,time,'linear','extrap'); % nT^2/Hz

B = B*10^-9; %  [T]
ne= ne*10^6; %  [m-3]
ni= ni*10^6; %  [m-3]

q = 1.602*10^-19; % [C]
me = 9.11*10^-31; % [kg]
epsilon = 8.85*10^-12; %[F.m-1]
mp = 1.67*10^-27; % [kg]
c = 3*10^8;  % [m/s]

%  Cyclotron frequency
omega_e = q*B/me; % rad/s

%  Plasma frequency
w_pe = (ne*q^2/(me*epsilon)).^0.5; % rad/s

% Certain Parameters
M = me/mp; 
Z = [1, 1, 1, 2]; % Positive charge on H+, O+, He+, He2+  
Beta = [1, 16, 4, 4];
fraction = [1.05,0.02,0.005,0.05]/sum([1.05,0.02,0.005,0.05]);
%  Using the paper: Lennartsson 1986, http://iopscience.iop.org/article/10.1088/0031-8949/36/2/029/pdf 
eta = repmat(fraction',1,1000).*repmat((ni./ne),4,1); % eta_i = n_i/n_e

% Calculating maximum cut off frequency
amp=1*10^0; %Power of the maximum cut off in nT^2/Hz
for i=1:1:length(time)
 [M1,I] =max(((wvpow(i,:)-amp)>0).*xwvfgm.freq);
 fuc(i)=xwvfgm.freq(I);
end;
w_uc = fuc*2*pi; % rad/s :  Maximum cut-off frequency

[M2,J]=max(w_uc);
[M2,J0]=min(w_uc);
% 
% for j=1:1:4
%     omega_i(j,:) = Z(j)*(omega_e*me/(Beta(j)*mp));
%     A(j,:) = (Beta(j)./(M*(Z(j)^2)*eta(j,:))).*((omega_e./w_pe).^2).*(1-w_uc./omega_i(j,:)).*(omega_i(j,:)./w_uc).^2;
%     A(A<0)=0;
%     E_min(j,:)=((A(j,:)+1).^0.5)*(me*c^2)*6.2415*10^18; %eV
% end;
I=168;
for j=1:1:4
    omega_i(j,:) = Z(j)*(omega_e*me/(Beta(j)*mp));
    w_uc_m=(0:0.00001:1)*(omega_i(1,I));
    A1(j,:) = (Beta(j)./(M*(Z(j)^2)*(eta(j,I)))).*((omega_e(1,I)./w_pe(1,I)).^2).*(1-w_uc_m./(omega_i(j,I))).*((omega_i(j,I))./w_uc_m).^2;
    E1_min(j,:)=((A1(j,:)).^0.5)*(me*c^2)*6.2415*10^18; %eV Kinetic energy
end;

%% Plotting
hFig=figure(5);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 60; %in mm
p.pack({{panelSize} {panelSize}},{{80}});
p.marginleft=30;
p.marginright=10;
p.fontsize=12;
p.de.margin=4;
p(1).margintop=10;
p(2).margintop=8;
p.select('all');

p(1,1).select();
ylabel({'Minimum cyclotron resonance','energy with EMIC waves','[MeV]'});
set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1],'XLim',[0 1],...
    'XTickLabel','',...
    'YLim',[0.01 100],'YTick',[0.01 0.1 1 10 100],'YTickLabel',{'0.01', '0.1','1','10','100'},...
    'YScale','log','XScale','linear');
xlabel('$\omega_{UC}$/$\Omega_{H^{+}}$','Interpreter','Latex');
grid on;

p(2,1).select();
ylabel({'$Wave\,power\,|B_{x}|$','\,\,\,\,\,$[nT/\sqrt{Hz}]$'},'Interpreter','Latex','FontSize',12);
set(gca,'XTick',[0,0.2,0.4,0.6,0.8,1],'XLim',[0 1],...
    'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'},...
    'YLim',10.^[-4 4],'YTick',10.^[-4 -3 -2 -1 0 1 2 3 4],...
    'YScale','log','XScale','linear');
xlabel('$\omega_{thm}$/$\Omega_{H^{+}}$','Interpreter','Latex','FontSize',12);
grid on;
%%

p(1,1).select();
plot(w_uc_m(1,:)./(omega_i(1,I)),E1_min(1,:)'/10^6,'-r');
hold on;
plot(w_uc_m(1,:)./(omega_i(1,I)),E1_min(2,:)'/10^6,'-k');
hold on;
plot(w_uc_m(1,:)./(omega_i(1,I)),E1_min(3,:)'/10^6,'--m');
hold on;
plot(w_uc_m(1,:)./(omega_i(1,I)),E1_min(4,:)'/10^6,'-.m');
hold on;
legend('H^+','O^+','He^+','He^2^+');

% h=subplot(3,1,2);
% g=area(time,w_uc./omega_i(1,:));
% g.FaceColor=[0.75 0 0];
% DateNumBeg=time(1);
% DateNumEnd=time(end);
% dt=0.5;
% TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd)];
% set(h,'XTick',TTick);
% set(h,'XTickLabel',{datestr(TTick,'HH:MM')});
% set(h,'XLim',[DateNumBeg  DateNumEnd]);
% set(get(h,'XLabel'),'String','Universal Time [HH:MM]');
% set(get(h,'YLabel'),'String','\omega_U_C/\Omega_p');

p(2,1).select();
tin=J; 
plot(2*pi*xwvfgm.freq./mean(omega_i(1,tin)),mean(wvpow(tin,:),1),'-r','LineWidth',2);
hold on;
plot(2*pi*xwvfgm.freq./mean(omega_i(1,J0)),mean(wvpow(J0,:),1),'-r','LineWidth',1);
legend(['After onset ',datestr(time(tin),'HH:MM'),' UT'],['Before onset ',datestr(time(J0),'HH:MM'),' UT']);

