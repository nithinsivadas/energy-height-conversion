%% Time averaged Energy Spectra before and after the onset
%  The routine creates the time-averaged energy spectra before and after
%  the substorm onset - from space and ground. 
%  Written by Nithin Sivadas, 29th June 2016. 

clear all;

load /home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/energy_spectra_13_june_2016.mat
% THEMIS-D Pitch  Angle
load /home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/pitch_angles.mat


data1.ssteflux=cell2mat(sst.eflux);
data1.esaeflux=cell2mat(esa.eflux);
data1.ssttime = sst.etime;
data1.esatime = esa.etime;
data1.pa=cell2mat(sst.epa);

% Solid angle in steridian
domega = 2*pi*(1-cos(11.5*pi/180));

% E_parallel/Total_Energy : Loss cone energy flux ratio
A_high = sum(data1.ssteflux(:,[1,8])*domega,2)./sum(data1.ssteflux(:,:)*domega,2);
A_mid  = A_high;
A_low  = sum(data1.esaeflux(:,[1,8])*domega,2)./sum(data1.esaeflux(:,:)*domega,2);

A_h = interp1(data1.ssttime,A_high,data_thm.time);
A_m = A_h;
A_l = interp1(data1.esatime,A_low,data_thm.time);
A_l(75) = A_l(76);
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
thm_B = mean(data_thm.E(75:114,:).*[repmat(A_h(75:114,:),1,length(data_thm.ebin)-low_t_ebin+1),repmat(A_h(75:114,:),1,low_t_ebin-1)],1);
thm_A = mean(data_thm.E(146:166,:).*[repmat(A_h(146:166,:),1,length(data_thm.ebin)-low_t_ebin+1),repmat(A_h(146:166,:),1,low_t_ebin-1)],1);


%% PFISR

% Splitting into high energy >300 keV, mid 30-300 keV and low energy <30 keV
pfisr = data_pfisr.E;

% Before time: 10:01 AM to 11:13 AM;  After time: 11!:46 UT to 12:17 UT
pfisr_B = mean(pfisr(57:92,:),1);
pfisr_A = mean(pfisr(108:123,:),1);

figure; loglog(data_thm.ebin,thm_B,'b'); hold on; loglog(data_thm.ebin, thm_A,'r');
hold on;
loglog(data_pfisr.ebin,pfisr_B,'.-b'); hold on; loglog(data_pfisr.ebin, pfisr_A,'.-r');

