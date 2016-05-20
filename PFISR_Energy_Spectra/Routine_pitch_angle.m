
load pitch_angles.mat

% Normalized pitch angle flux
data.ssteflux=cell2mat(sst.eflux)./repmat(sum(cell2mat(sst.eflux)*(22.5/180),2),1,8);
data.sstiflux=cell2mat(sst.iflux)./repmat(sum(cell2mat(sst.iflux)*(22.5/180),2),1,8);
data.esaiflux=cell2mat(esa.iflux)./repmat(sum(cell2mat(esa.iflux)*(22.5/180),2),1,8);
data.esaeflux=cell2mat(esa.eflux)./repmat(sum(cell2mat(esa.eflux)*(22.5/180),2),1,8);
data.pa=cell2mat(sst.epa);

figure;

%  In the form of plot_eflux function

%% SST Electron
data_se.E=10.^(data.ssteflux);
data_se.ebin=data.pa;
data_se.time=sst.etime;

subplot(4,1,3)
plot_eflux(data_se,0.5,'SST e^- PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% SST Ion
data_si.E=10.^(data.sstiflux);
data_si.ebin=data.pa;
data_si.time=sst.itime;

subplot(4,1,1)
plot_eflux(data_si,0.5,'SST i^+ PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% ESA Electron
data_ee.E=(10.^data.esaeflux);
data_ee.ebin=data.pa;
data_ee.time=esa.etime;

subplot(4,1,4)
plot_eflux(data_ee,0.5,'ESA e^- PAD');
cmin=0;cmax=1.5;
set(gca,'YScale','linear');
ylabel('Pitch angle in Degrees');
caxis([cmin cmax]);

%% ESA Ion
data_ei.E=10.^(data.esaiflux);
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
subplot(4,1,1)
plot(data_si.time,sum(data_si.E(:,[3,4,5,6]),2)./sum(data_si.E(:,[1,2,7,8]),2));
dt=0.5;
DateNumBeg=min(data_si.time);
DateNumEnd=max(data_si.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','E_/perp/E_/parallel');
ylim([0 2]);
title('SST Ion Anisotropy');

subplot(4,1,2)
plot(data_ei.time,sum(data_ei.E(:,[3,4,5,6]),2)./sum(data_ei.E(:,[1,2,7,8]),2));
dt=0.5;
DateNumBeg=min(data_ei.time);
DateNumEnd=max(data_ei.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','E_/perp/E_/parallel');
ylim([0 2]);
title('ESA Ion Anisotropy');

subplot(4,1,3)
plot(data_se.time,sum(data_se.E(:,[3,4,5,6]),2)./sum(data_se.E(:,[1,2,7,8]),2));
dt=0.5;
DateNumBeg=min(data_se.time);
DateNumEnd=max(data_se.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','E_/perp/E_/parallel');
ylim([0 2]);
title('SST Electron Anisotropy');

subplot(4,1,4)
plot(data_ee.time,sum(data_ee.E(:,[3,4,5,6]),2)./sum(data_ee.E(:,[1,2,7,8]),2));
dt=0.5;
DateNumBeg=min(data_ee.time);
DateNumEnd=max(data_ee.time);
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd),DateNumEnd];
% set(gca, 'Layer','top') 
set(gca,'XTick',TTick);
set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','E_/perp/E_/parallel');
ylim([0 2]);
title('ESA Electron Anisotropy');
