%% Pitch angle anisotropy
%  Plot the pitch angle distribution of energies  <30 keV and <30 keV of 
%  electrons and ions using THEMIS-D. 

load '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/pitch_angle_8_to_13_UT.mat'

% Normalized pitch angle flux
data.ssteflux=(sst.eflux)./repmat(sum((sst.eflux)*(22.5/180),2),1,8);
data.sstiflux=(sst.iflux)./repmat(sum((sst.iflux)*(22.5/180),2),1,8);
data.esaiflux=(esa.iflux)./repmat(sum((esa.iflux)*(22.5/180),2),1,8);
data.esaeflux=(esa.eflux)./repmat(sum((esa.eflux)*(22.5/180),2),1,8);
data.pa=cell2mat(sst.epa);
% 
% data.ssteflux=cell2mat(sst.eflux);
% data.sstiflux=cell2mat(sst.iflux); 
% data.esaiflux=cell2mat(esa.iflux);
% data.esaeflux=cell2mat(esa.eflux);


figure;

%  In the form of plot_eflux function

%% SST Electron
data_se.E=(data.ssteflux);
data_se.ebin=data.pa;
data_se.time=sst.etime;

subplot(4,1,3)
plot_eflux(data_se,0.5,'SST e^- PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% SST Ion
data_si.E=(data.sstiflux);
data_si.ebin=data.pa;
data_si.time=sst.itime;

subplot(4,1,1)
plot_eflux(data_si,0.5,'SST i^+ PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% ESA Electron
data_ee.E=(data.esaeflux);
data_ee.ebin=data.pa;
data_ee.time=esa.etime;

subplot(4,1,4)
plot_eflux(data_ee,0.5,'ESA e^- PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% ESA Ion
data_ei.E=(data.esaiflux);
data_ei.ebin=data.pa;
data_ei.time=esa.itime;

subplot(4,1,2)
plot_eflux(data_ei,0.5,'ESA i^+ PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% SST Ion Anisotropy
figure;
subplot(2,1,1)
plot(data_si.time,sum(data_si.E(:,[3,4,5,6]),2)./sum(data_si.E(:,[1,2,7,8]),2)-1,'r');
dt=0.5;
DateNumBeg=min(data_si.time);
DateNumEnd=max(data_si.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','(\phi_\perp / \phi_|_|) - 1');
ylim([-1 1]);
% label();
title('Ion Anisotropy');
% subplot(4,1,2)
hold on;
plot(data_ei.time,sum(data_ei.E(:,[3,4,5,6]),2)./sum(data_ei.E(:,[1,2,7,8]),2)-1);
dt=0.5;
% DateNumBeg=min(data_ei.time);
% DateNumEnd=max(data_ei.time);
% TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
% set(gca,'XTick',TTick);
% set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
% set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
% set(get(gca,'YLabel'),'String','(\phi_\perp / \phi_|_|) - 1');
% ylim([-1 1]);
legend('Ions > 25 keV','Ions < 25 keV');

subplot(2,1,2)
plot(data_se.time,sum(data_se.E(:,[3,4,5,6]),2)./sum(data_se.E(:,[1,2,7,8]),2)-1,'r');
dt=0.5;
DateNumBeg=min(data_se.time);
DateNumEnd=max(data_se.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','(\phi_\perp / \phi_|_|) - 1');
ylim([-1 1]);
title('Electron Anisotropy');

hold on;
% subplot(4,1,4)
plot(data_ee.time,sum(data_ee.E(:,[3,4,5,6]),2)./sum(data_ee.E(:,[1,2,7,8]),2)-1);
% dt=0.5;
% DateNumBeg=min(data_ee.time);
% DateNumEnd=max(data_ee.time);
% TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% % set(gca, 'Layer','top') 
% set(gca,'XTick',TTick);
% set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
% set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
% set(get(gca,'YLabel'),'String','(\phi_\perp / \phi_|_|) - 1');
% ylim([-1 1]);
% title('Electron Anisotropy < 30 keV');
legend('Electrons > 30 keV', 'Electrons < 30 keV');
