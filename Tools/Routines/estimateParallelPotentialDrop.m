%% Estimating Parallel Potential drop

clear all;
data=get_all_time_series_data();
load('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')

pad_low = thd.pitchAngleDistribution.esa.flux;
pa =cell2mat(thd.pitchAngleDistribution.esa.pitchAngle);
pa_time = thd.pitchAngleDistribution.esa.time;
pad_high = thd.pitchAngleDistribution.sst.flux;

loss_cone_angle = 1.5; %deg
%%
for thisTime=1:1:length(pa_time)
  lc_flux_low(thisTime)=interp1(pa,pad_low(thisTime,:),loss_cone_angle,'linear','extrap');
  avg_flux_low(thisTime)=mean(pad_low(thisTime,:));
  ratio_flux_low(thisTime)=lc_flux_low(thisTime)./avg_flux_low(thisTime);
  
  lc_flux_high(thisTime)=interp1(pa,pad_high(thisTime,:),loss_cone_angle,'linear','extrap');
  avg_flux_high(thisTime)=mean(pad_high(thisTime,:));
  ratio_flux_high(thisTime)=lc_flux_high(thisTime)./avg_flux_high(thisTime);
end;

%%
thd_eflux=data(8).zValue;
thd_eflux(1:32,:)=thd_eflux(1:32,:).*repmat(ratio_flux_low,32,1);
thd_eflux(33:end,:)=thd_eflux(33:end,:).*repmat(ratio_flux_high,10,1);
thd_time = data(8).time;
thd_eBin = data(8).yAxis;

pfisr_eflux = data(5).zValue;
pfisr_time = data(5).time;
pfisr_eBin = data(5).yAxis;

%%
for thisTime=1:1:length(thd_time)
     [c,I]=max(thd_eflux(:,thisTime));
     sample_Y=log(thd_eflux(:,thisTime));
     sample_X = thd_eBin';
     outliers=excludedata(sample_X,sample_Y,'box',[thd_eBin(I) thd_eBin(end) log(10.^8) log(10.^13)]);
     [fit,gof]= fit(sample_X,sample_Y,'poly1','Exclude',outliers);
     a_T(thisTime)=fit.p1;
     b_T(thisTime)=fit.p2;
     rsquare_T(thisTime)=gof.rsquare;
     peak_energy_T(thisTime)=thd_eBin(I);
     clearvars fit gof
     
     timeNo=find_time(pfisr_time, datestr(thd_time(thisTime))); 
     [c,I]=max(pfisr_eflux(:,timeNo));
     sample_Y1=log(pfisr_eflux(:,timeNo));
     sample_X1= pfisr_eBin;
     outliers1=excludedata(sample_X1,sample_Y1,'box',[pfisr_eBin(I)' pfisr_eBin(end)' log(10.^8) log(10.^13)]);
     [fit1,gof1]= fit(sample_X1,sample_Y1,'poly1','Exclude',outliers1);
     a_P(thisTime)=fit1.p1;
     b_P(thisTime)=fit1.p2;
     rsquare_P(thisTime)=gof1.rsquare;
     peak_energy_P(thisTime)=pfisr_eBin(I);  
     
end;

% E0 = (b_P-b_T).*(-1*a_P.^-1); % Parallel Potential Drop
E0 = peak_energy_P-peak_energy_T; % Parallel Potential Drop
% figure; plot(thd_time,peak_energy_P/1000); hold on; plot(thd_time,peak_energy_T/1000); hold on;plot(thd_time,(peak_energy_P-peak_energy_T)/1000); label_time_axis(thd_time,true);
parallel_potential.V = E0/1000;
parallel_potential.time = thd_time;
parallel_potential.units = '[kV]';

