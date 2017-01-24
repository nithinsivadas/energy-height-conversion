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
     sample_Y=real(log10(thd_eflux(:,thisTime)));
     sample_X = thd_eBin';
     peak_energy_T(thisTime)=thd_eBin(I);
     [pks, locs, w, p] = findpeaks(sample_Y, sample_X);
     thisPeak = 1;
     while thisPeak<=length(pks)
         if(pks(thisPeak)>=pks(1)-1)
            local.pk(1,thisTime) = pks(thisPeak);
            local.loc(1,thisTime) = locs(thisPeak);
            local.pkno(1,thisTime) = thisPeak;
         end
         thisPeak = thisPeak+1;
     end
        
     timeNo=find_time(pfisr_time, datestr(thd_time(thisTime))); 
     [c1,I1]=max(pfisr_eflux(:,timeNo));
     sample_Y1=log10(pfisr_eflux(:,timeNo));
     sample_X1= pfisr_eBin;
     peak_energy_P(thisTime)=pfisr_eBin(I1);
     clearvars pks locs w p;
     [pks, locs, w, p] = findpeaks(sample_Y1, sample_X1);
     thisPeak = 1;
     while thisPeak<=length(pks)
         if(pks(thisPeak)>=pks(1)-1)
            local.pk(2,thisTime) = pks(thisPeak);
            local.loc(2,thisTime) = locs(thisPeak);
            local.pkno(2,thisTime) = thisPeak;
         end
         thisPeak = thisPeak+1;
     end
     
end;

% E0 = (b_P-b_T).*(-1*a_P.^-1); % Parallel Potential Drop
% E0 = peak_energy_P-peak_energy_T; % Parallel Potential Drop
E0 = local.loc(2,:)-local.loc(1,:);
parallel_potential.V = E0/1000;
parallel_potential.time = thd_time;
parallel_potential.units = '[kV]';
parallel_potential.peakinfo=local;

% figure; plot(thd_time,local.loc(2,:)/1000); hold on; plot(thd_time,local.loc(1,:)/1000); hold on;
figure; plot(thd_time,(E0)/1000); label_time_axis(thd_time,true);
ylabel('Potential Drop [kV]');


