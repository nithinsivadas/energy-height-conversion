clear all;

load /home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/energy_spectra_13_june_2016.mat


n=3;

% THEMIS
low_t_ebin=32; %30 keV
mid_t_ebin=36+n; %300 keV

% PFISR
low_p_ebin=32; %30 keV
mid_p_ebin=36+n; %300 keV


%% THEMIS
thm = data_thm.E;

% Splitting in to before onset and after onset Befgore: 10:00 - 11:02;
% After : 11:54-12:25
thm_B = mean(data_thm.E(75:114,:),1);
thm_A = mean(data_thm.E(146:166,:),1);


%% PFISR


% Splitting into high energy >300 keV, mid 30-300 keV and low energy <30 keV
pfisr = data_pfisr.E;

% Before time: 10:01 AM to 11:13 AM;  After time: 11!:46 UT to 12:17 UT
pfisr_B = mean(pfisr(57:92,:),1);
pfisr_A = mean(pfisr(108:123,:),1);

figure; loglog(data_thm.ebin,thm_B,'b'); hold on; loglog(data_thm.ebin, thm_A,'r');
hold on;
loglog(data_pfisr.ebin,pfisr_B,'.-b'); hold on; loglog(data_pfisr.ebin, pfisr_A,'.-r');

