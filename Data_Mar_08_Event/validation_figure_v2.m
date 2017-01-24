% Plot production rates 
% 1. estimated by PFISR
% 2. estimated by the forward model using the energy spectra obtained from
% inversion
% 3. estimated by the forward model using the energy spectra obtained from
% THEMIS
function [] = validation_figure_v2( date_string )

% load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/raw_pfisr_thm_production_rates_v2.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/A_matrix_production_rates.mat');
load('/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/sample_output_20_jul_2016.mat');
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
A =generate_A(data.h,data_thm.ebin(11:end)');
% A=A_thm_forwd(:,:);

lc_eflux(:,1:31) = data_thm.E(:,1:31).*repmat(lc_ratio_sst,1,31);
lc_eflux(:,32:42) = data_thm.E(:,32:42).*repmat(lc_ratio_sst,1,11);

% q_forwd_thm=A*(2*pi*data_thm.E(i_TE,9:end)'./data_thm.ebin(9:end)');
q_forwd_thm=A*(2*pi*lc_eflux(i_TE,11:end)'./data_thm.ebin(11:end)');

subplot(1,2,1)
plot(data.q_PFISR(:,i_pq),data.h);
% hold on;
% plot(data.q_PFISR(:,i_pq)+data.dn_dt(:,i_pq),data.h,'g');
hold on; 
plot(data.q_THEMIS(:,i_pE),data.h,'--r');
hold on;
plot(q_forwd_thm,data.h,'-.b');

set(gca,'XScale','log','YLim',[60 200]); 
legend(['Q_P_F_I_S_R'],...
    ['Q_f_i_t=A\phi_e_s_t'],...
    ['Q_t_h_m=A\phi_o_b_s'],'location','northwest');
xlabel('Production rate [m^-^3 s^-^1]');
ylabel('Altitude [km]');

subplot(1,2,2)
plot(data_pfisr.ebin/1000,data_pfisr.E(i_pE,:),'--r');
% hold on;
% plot(data_pfisr_new.E(i_pE,:),data_pfisr.ebin);
hold on; 
plot(data_thm.ebin/1000,data_thm.E(i_TE,:),'-.b');
set(gca,'XScale','log','XLim',[10^0 10^3],'YScale','log'); 
legend(['\phi_f_i_t'],...
      ['\phi_t_h_m'],'location','northeast');
xlabel('Energy [keV]');
ylabel('Differential Energy Flux [eV/cm^2  sr s eV]'); 
title(datestr(data.time_q(i_pq)));

end