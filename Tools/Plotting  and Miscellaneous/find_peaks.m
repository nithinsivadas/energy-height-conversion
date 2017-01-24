function [ local, sample_Y, sample_X, sample_Y1, sample_X1 ] = find_peaks( thisTimeStr, thd_time, thd_eflux, thd_eBin, pfisr_time, pfisr_eflux, pfisr_eBin )
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

thisTime=find_time(thd_time, thisTimeStr);

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
     
     figure;
     plot(sample_X,sample_Y,'k.-');
     hold on;
     plot(locs,pks,'vk');
     hold on;
     plot(local.loc(1,thisTime),local.pk(1,thisTime)+0.2,'*k');
     hold on;
     
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
 
     plot(sample_X1, sample_Y1, 'r.-');
     hold on;
     plot(locs,pks,'vr');
     hold on;
     plot(local.loc(2,thisTime),local.pk(2,thisTime)+0.2,'*r');
          
     ylim([8 13]);
     ylabel('log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
     xlabel('Energy [eV]');
     set(gca,'xScale','log');
     title([datestr(thd_time(thisTime),'HH:MM'),' UT ',num2str((local.loc(2,thisTime)-local.loc(1,thisTime))/1000),' kV']);
     legend('THEMIS Energy Spectra','THM Peaks','THM Chosen Peak','PFISR Energy Spectra','PFISR Peaks','PFISR Chosen Peak','location','eastoutside');

end

