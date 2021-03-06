%% Energy Spectra measured by THEMIS and PFISR
% May 17th, 2016

% Required data files
%
%

clear all;

DateNumBeg=datenum('26-Mar-2008 08:00');
DateNumEnd=datenum('26-Mar-2008 13:00');

%% Acquiring THEMIS Energy Spectra

data_esa = read_tplot_ascii('Data/thd_peef_en_eflux_20080326.txt','Data/thd_peef_en_eflux_v_20080326.txt');
data_sst = read_tplot_ascii('Data/thd_psef_en_eflux_20080326.txt','Data/thd_psef_en_eflux_v_20080326.txt');
data_thm = combine_energies(data_esa,data_sst);

% data_thm -> E, ebin, time

%% Acquiring PFSIR Energy Spectra

[data]=converth2Ev4(data_thm, 'Data/DataFile_2008_1.h5', DateNumBeg, DateNumEnd); % h values is the default value from PFISR data DataFile_2008_1.h5

% Converting the data from PFISR into plottable format
data_pfisr.E=data.eflux';
data_pfisr.ebin=data.ebin;
data_pfisr.time=data.time_PFISR_E;
data_thm.E(data_thm.E==0)=10^1;
figure;
subplot(2,1,1)
plot_eflux(data_pfisr,0.5,'Energy Spectra from PFISR 26 Mar 2008');
subplot(2,1,2)
plot_eflux(data_thm,0.5,'Energy Spectra from THEMIS-D 26 Mar 2008');

save('Data/sample_output_20_jul_2016.mat');
