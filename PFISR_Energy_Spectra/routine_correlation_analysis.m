%% Correlation  Analysis
clear all;

load /home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/sample_output_19_may_2016.mat
% T1=4;
T1=145;
% Interpolating
data_pfisr_i.E=interp2(data_pfisr.time,data_pfisr.ebin,data_pfisr.E',data_thm.time,data_thm.ebin(9:end))';
data_pfisr_i.E=data_pfisr_i.E(T1:end-3,:);
data_pfisr_i.ebin=data_thm.ebin(9:end);
data_pfisr_i.time=data_thm.time(T1:end-3);

data_thm_i.E=data_thm.E(T1:end-3,9:end);
data_thm_i.ebin=data_thm.ebin(9:end);
data_thm_i.time=data_thm.time(T1:end-3);

for i=1:1:34
    x= data_thm_i.E(:,i);
    z= data_pfisr_i.E(:,i);
%     z(1:135)=0;  %11:40:59 UT
%     z(145:end)=0;%11:57:01 UT
%     [r(:,i), lags(i,:)] = xcorr(x,z,'coeff');
    [R1,P1] = corrcoef(x,z);                                                                                                                                                
    r(i) = R1(1,2);
    p(i) = P1(1,2);
end;

% figure; plot(lags(1,:),r(:,:));
figure;plot(data_thm_i.ebin/1000,r);

figure;
subplot(2,1,1)
plot_eflux(data_pfisr_i);
subplot(2,1,2)
plot_eflux(data_thm_i);

% Correlation

