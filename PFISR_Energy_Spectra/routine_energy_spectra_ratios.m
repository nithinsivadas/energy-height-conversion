%% Ratio of energy spectra before and after onset

clear all;

load('/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/sample_output_20_jul_2016.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/thd_pa.mat');

lc_ratio_esa=sum(thd_pa.esa.flux(:,[1,8]),2)./sum(thd_pa.esa.flux,2);
lc_ratio_sst=sum(thd_pa.sst.flux(:,[1,8]),2)./sum(thd_pa.sst.flux,2);

time.thm.bo(1)=datenum('2008-03-26/10:00:00');
time.thm.bo(2)=datenum('2008-03-26/11:02:00');
time.thm.ao(1)=datenum('2008-03-26/11:54:00');
time.thm.ao(2)=datenum('2008-03-26/12:25:00');

time.pfisr.bo(1)=datenum('2008-03-26/10:00:00');
time.pfisr.bo(2)=datenum('2008-03-26/11:02:00');
time.pfisr.ao(1)=datenum('2008-03-26/11:46:00');
time.pfisr.ao(2)=datenum('2008-03-26/12:17:00');

[M,time.thm.i_bo(1)] = min(abs(data_thm.time-time.thm.bo(1))); 
[M,time.thm.i_bo(2)] = min(abs(data_thm.time-time.thm.bo(2)));
[M,time.thm.i_ao(1)] = min(abs(data_thm.time-time.thm.ao(1))); 
[M,time.thm.i_ao(2)] = min(abs(data_thm.time-time.thm.ao(2)));

[M,time.pfisr.i_bo(1)] = min(abs(data_pfisr.time-time.pfisr.bo(1))); 
[M,time.pfisr.i_bo(2)] = min(abs(data_pfisr.time-time.pfisr.bo(2)));
[M,time.pfisr.i_ao(1)] = min(abs(data_pfisr.time-time.pfisr.ao(1))); 
[M,time.pfisr.i_ao(2)] = min(abs(data_pfisr.time-time.pfisr.ao(2)));

eflux.thm.bo=mean(data_thm.E(time.thm.i_bo(1):1:time.thm.i_bo(2),:).*repmat(lc_ratio_sst(time.thm.i_bo(1):1:time.thm.i_bo(2)),1,length(data_thm.ebin)),1);
eflux.thm.ao=mean(data_thm.E(time.thm.i_ao(1):1:time.thm.i_ao(2),:).*repmat(lc_ratio_sst(time.thm.i_ao(1):1:time.thm.i_ao(2)),1,length(data_thm.ebin)),1);

eflux.pfisr.bo=mean(data_pfisr.E(time.pfisr.i_bo(1):1:time.pfisr.i_bo(2),:),1);
eflux.pfisr.ao=mean(data_pfisr.E(time.pfisr.i_ao(1):1:time.pfisr.i_ao(2),:),1);

figure;
subplot(2,1,1)
plot(data_thm.ebin,eflux.thm.bo,'-b'); hold on; plot(data_thm.ebin, eflux.thm.ao,'--b');
set(gca,'xscale','log');
set(gca,'yscale','log');
subplot(2,1,2)
plot(data_pfisr.ebin,eflux.pfisr.bo,'-r'); hold on; plot(data_pfisr.ebin, eflux.pfisr.ao,'--r');
set(gca,'xscale','log');
set(gca,'yscale','log');
