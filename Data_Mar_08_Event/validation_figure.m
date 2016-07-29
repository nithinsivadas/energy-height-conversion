% Plot production rates 
% 1. estimated by PFISR
% 2. estimated by the forward model using the energy spectra obtained from
% inversion
% 3. estimated by the forward model using the energy spectra obtained from
% THEMIS
function [] = validation_figure( date_string )

load('/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/energy_spectra_13_june_2016.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/raw_pfisr_thm_production_rates.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/thd_pa.mat');
figure;

lc_ratio_esa=sum(thd_pa.esa.flux(:,[1,8]),2)./sum(thd_pa.esa.flux,2);
lc_ratio_sst=sum(thd_pa.sst.flux(:,[1,8]),2)./sum(thd_pa.sst.flux,2);

select_time=datenum(date_string);

[M,i_pq] = min(abs(data.time_q-select_time)); %Time index for PFISR production rate data
[M,i_pE] = min(abs(data.time_PFISR_E-select_time)); %Time index for PFISR derived energy spectral data
[M,i_TE] = min(abs(data.time_E-select_time)); %Time index for THEMIS energy spectral data
[M,i_pn] = min(abs(data.time_ne-select_time)); %Time index for PFISR electron density data

% altitude=data.h;
% A =generate_A(data.h,data_thm.ebin');
A=A_thm_forwd(:,9:end);

lc_eflux(:,1:31) = data_thm.E(:,1:31).*repmat(lc_ratio_sst,1,31);
lc_eflux(:,32:42) = data_thm.E(:,32:42).*repmat(lc_ratio_sst,1,11);

% q_forwd_thm=A*(2*pi*data_thm.E(i_TE,9:end)'./data_thm.ebin(9:end)');
q_forwd_thm=A*(2*pi*lc_eflux(i_TE,9:end)'./data_thm.ebin(9:end)');

plot(data.q_PFISR(:,i_pq),data.h); 
hold on; 
plot(data.q_THEMIS(:,i_pE),data.h,'--r');
hold on;
plot(q_forwd_thm,data.h,'-.b');

set(gca,'XScale','log','YLim',[60 250]); 
legend([datestr(data.time_q(i_pq),'HH:MM:SS'),' UT Q_P_F_I_S_R'],...
    [datestr(data.time_PFISR_E(i_pE),'HH:MM:SS'),' UT Q_n_e_w=A\phi_e_s_t'],...
    [datestr(data.time_E(i_TE),'HH:MM:SS'),' UT Q_t_h_m=A\phi_o_b_s'],'location','northwest');

end