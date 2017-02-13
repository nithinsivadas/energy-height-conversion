function [ local, thdEflux, thdeBin, pfisrEflux, pfisreBin ] = find_peaks( thisTimeStr, thd_time, thd_eflux, thd_eBin, pfisr_time, pfisr_eflux, pfisr_eBin )
%% find_peaks.m Identify peaks in the energy spectra of THEMIS and PFISR
%--------------------------------------------------------------------------
% Input
%------
% thisTimeStr - String of the time instant of interest e.g. 26 Mar 2008
%               11:00 am
% thd_time    - Time array in matlab units for THEMIS data [nTx1}
% thd_eflux   - 2-D energy flux matrix [nExnT]
% thd_eBin    - Energy bin array [nEx1]
% pfisr_time  - Time array in matlab units for PFISR data [nTx1]
% pfisr_eflux - 2-D energy flux matrix [nExnT]
% pfisr_eBin  - Energy bin array [nEx1]
%--------------------------------------------------------------------------
% Output
%------
% local       - Structure containg details of the local peaks in the energy spectra 
%      ->pk   - Energy flux value of the peak
%      ->locs - Location of the peak / energy value
%      ->pkno - Number of the peak
% thdEflux    - An array of themis eflux values at thisTime
% thdeBin     - An array of themis energy bin values at thisTime
% pfisrEflux    - An array of PFISR eflux values at thisTime
% pfisreBin     - An array of PFISR energy bin values at thisTime
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
%% THEMIS
    thisTime=find_time(thd_time, thisTimeStr);

     [c,I]=max(thd_eflux(:,thisTime));
     thdEflux=real(log10(thd_eflux(:,thisTime)));
     thdeBin = thd_eBin';
     peak_energy_T(thisTime)=thd_eBin(I);
     [pks, locs, w, p] = findpeaks(thdEflux, thdeBin);
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
     plot(thdeBin,thdEflux,'k.-');
     hold on;
     plot(locs,pks,'vk');
     hold on;
     plot(local.loc(1,thisTime),local.pk(1,thisTime)+0.2,'*k');
     hold on;

%% PFISR
     timeNo=find_time(pfisr_time, datestr(thd_time(thisTime))); 
     [c1,I1]=max(pfisr_eflux(:,timeNo));
     pfisrEflux=log10(pfisr_eflux(:,timeNo));
     pfisreBin= pfisr_eBin;
     peak_energy_P(thisTime)=pfisr_eBin(I1);
     clearvars pks locs w p;
     [pks, locs, w, p] = findpeaks(pfisrEflux, pfisreBin);
     thisPeak = 1;
     while thisPeak<=length(pks)
         if(pks(thisPeak)>=pks(1)-1)
            local.pk(2,thisTime) = pks(thisPeak);
            local.loc(2,thisTime) = locs(thisPeak);
            local.pkno(2,thisTime) = thisPeak;
         end
         thisPeak = thisPeak+1;
     end
 
     plot(pfisreBin, pfisrEflux, 'r.-');
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

