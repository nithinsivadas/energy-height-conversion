%% Ratio of time averaged energy flux
%  Ratio of time-averaged energy flux before and after the onset

clear all;

load Data/energy_spectra_13_june_2016.mat
load Data/pitch_angles.mat
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

A_h = ones(size(A_h));
A_m = A_h;
A_l = ones(size(A_h));

n=0;
% THEMIS
% low_t_ebin=32; %30 keV
low_t_ebin=28; %10 keV
mid_t_ebin=36+n; %300 keV

% PFISR
% low_p_ebin=17; %30 keV
low_p_ebin=14; %10 keV
mid_p_ebin=20+n; %300 keV


%% THEMIS

% Splitting into high energy >300 keV, med 30-300 keV and low energy <30 keV
thm_high=sum(data_thm.E(:,mid_t_ebin:end).*repmat(A_h,1,length(data_thm.ebin)-mid_t_ebin+1).*repmat(diff(data_thm.ebin(1,(mid_t_ebin-1):end)),187,1),2);
thm_mid =sum(data_thm.E(:,low_t_ebin:mid_t_ebin-1).*repmat(A_m,1,mid_t_ebin-low_t_ebin).*repmat(diff(data_thm.ebin(1,low_t_ebin:mid_t_ebin)),187,1),2);
thm_low =sum(data_thm.E(:,1:low_t_ebin-1).*repmat(A_l,1,low_t_ebin-1).*repmat(diff(data_thm.ebin(1,1:low_t_ebin)),187,1),2);

% Splitting in to before onset and after onset Befgore: 10:00 - 11:02;
% After : 11:54-12:25
thm_high_B = mean(thm_high(75:114));
thm_high_A = mean(thm_high(146:166));
ratio_thm_high=thm_high_A/thm_high_B

thm_mid_B = mean(thm_mid(75:114));
thm_mid_A = mean(thm_mid(146:166));
ratio_thm_mid=thm_mid_A/thm_mid_B

thm_low_B = mean(thm_low(76:114));
thm_low_A = mean(thm_low(146:166));
ratio_thm_low=thm_low_A/thm_low_B

%% PFISR


% Splitting into high energy >300 keV, mid 30-300 keV and low energy <30 keV
pfisr_high=sum(data_pfisr.E(:,mid_p_ebin:end).*repmat(diff(data_pfisr.ebin(1,mid_p_ebin-1:end)),141,1),2);
pfisr_mid =sum(data_pfisr.E(:,low_p_ebin:mid_p_ebin-1).*repmat(diff(data_pfisr.ebin(1,low_p_ebin:mid_p_ebin)),141,1),2);
pfisr_low =sum(data_pfisr.E(:,1:low_p_ebin-1).*repmat(diff(data_pfisr.ebin(1,1:low_p_ebin)),141,1),2);

% Before time: 10:01 AM to 11:13 AM;  After time: 11!:46 UT to 12:17 UT
pfisr_high_B = mean(pfisr_high(57:92));
pfisr_high_A = mean(pfisr_high(108:123));
ratio_pfisr_high=pfisr_high_A/pfisr_high_B

pfisr_mid_B = mean(pfisr_mid(57:92));
pfisr_mid_A = mean(pfisr_mid(108:123));
ratio_pfisr_mid=pfisr_mid_A/pfisr_mid_B

pfisr_low_B = mean(pfisr_low(57:92));
pfisr_low_A = mean(pfisr_low(108:123));
ratio_pfisr_low=pfisr_low_A/pfisr_low_B

y_pfisr =[ratio_pfisr_low, ratio_pfisr_mid, ratio_pfisr_high];
y_thm   =[ratio_thm_low, ratio_thm_mid, ratio_thm_high];
x = ['<30 keV', '30 - 300 keV', '300 keV'];

figure1 = figure; 
axes1 = axes('Parent',figure1,'YScale','log','XTick',[1 2 3]);

box(axes1,'on');
hold(axes1,'on');

bar1=bar([y_pfisr',y_thm']); 
set(bar1(1),'BaseValue',1);
title([num2str(data_pfisr.ebin(1,low_p_ebin)/1000),'-',num2str(data_pfisr.ebin(1,mid_p_ebin)/1000),'keV -PFISR ',num2str(data_thm.ebin(1,low_t_ebin)/1000),'-',num2str(data_thm.ebin(1,mid_t_ebin)/1000),'keV - THEMIS']);

display('PFISR Mid Bin');
data_pfisr.ebin(1,mid_p_ebin)/1000
display('THEMIS Mid Bin');
data_thm.ebin(1,mid_t_ebin)/1000
 



