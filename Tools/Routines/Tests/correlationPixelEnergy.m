clear all;
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-energyFlux.h5';
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20101018.001_bc_2min-energyFlux.h5';
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20080326.001_bc_15sec-energyFlux_v85.h5';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';

dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');
time = dascData.time;

timeMinIndx = find_time(time,'26-Mar-2008 11:00');
timeMaxIndx = find_time(time,'26-Mar-2008 11:45');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');

ibeam = 13;
energySlice = 1:1:500; %keV

%%
k=1;
probeLat = pfisrData.latitude(ibeam,1);
probeLon = pfisrData.longitude(ibeam,1);
multiWaitbar('Calculating...',0);
id = 1/(timeMaxIndx-timeMinIndx+1);
for itime = timeMinIndx:timeMaxIndx
    try
    tempData = get_2D_plot_inputs_at_time(h5FileStr,...
        'plotModeStr','OpticalImage',...
        'thisTimeIndx',itime);
    image=tempData.image;
    latitude = dascData.latitude;
    longitude = dascData.longitude;
    nanFlags = isnan(latitude) | isnan(longitude) | isnan(image);
    latitude(nanFlags) = [];
    longitude(nanFlags) = [];
    image(nanFlags) = [];
    F= scatteredInterpolant(latitude(:),longitude(:),image(:),'linear','none');
%     tempPixel = [F(probeLat,probeLon), F(probeLat,probeLon+0.05), F(probeLat+0.05,probeLon),...
%         F(probeLat-0.05,probeLon), F(probeLat,probeLon-0.5), F(probeLat,probeLon)];
%     pixels(k)=sum(tempPixel(:))./5;
    pixels(k) = F(probeLat,probeLon);
    thisTimeDASC(k)=tempData.thisTime;
    
    timePfisrIndx = find_time(pfisrData.time,datestr(dascData.time(itime))); 
    tempPfisrData = get_2D_plot_inputs_at_time(h5FileStr,...
        'plotModeStr','EnergyFluxMap',...
        'thisTimeIndx',timePfisrIndx);
    energyFlux(:,k) = interp1(pfisrData.zEnergyBin(ibeam,:),tempPfisrData.diffEnergyFlux(ibeam,:),energySlice(:)*1000);
%     energyFlux(k) = trapz(pfisrData.zEnergyBin(ibeam,:),tempPfisrData.diffEnergyFlux(ibeam,:));
    thisTimePFISR(k)=tempPfisrData.thisTime;   
    multiWaitbar('Calculating...','Increment',id);    
    k=k+1;
    catch ME;
    end   
end
multiWaitbar('Calculating...',1);    
multiWaitbar('CLOSEALL');

%%
timeMinIndx = find_time(thisTimeDASC,'26-Mar-2008 11:17');
timeMaxIndx = find_time(thisTimeDASC,'26-Mar-2008 11:22');

% for iTime = 1:1:100
    iTime=1;
    iTimeIndx = timeMinIndx+iTime-1:1:timeMaxIndx+iTime-1;
    thisTime = thisTimeDASC(iTimeIndx);
    % Time-corrected pixels and energyFluxes
    pixels1 = pixels(iTimeIndx);
    pixels1(pixels1>450)= nan;
    pixels1(pixels1<300)= nan;
    energyFlux(isnan(energyFlux)) = 10^6;
    timeNan(iTime) = any(isnan(pixels1));
    energyFlux1 = interp1(thisTimePFISR,energyFlux',thisTime,'linear','extrap')';
    energyFluxLog1 = log10(energyFlux1);

    % Normalize both
    energyFluxLog2 = normalize_matrix(energyFluxLog1);
    pixels1(isnan(pixels1))=nanmean(pixels1);
    pixels1 = log10(pixels1);
    pixels2 = normalize_matrix(pixels1);

    % Find correlation
    r = zeros(1,length(energySlice));
    for i=1:1:length(energySlice)
        r(i)=xcorr(pixels2,energyFluxLog2(i,:),0,'coeff');
    end
%     
%     [~,indx] = max(r);
%     peakEnergy(iTime) = energySlice(indx);
%     peakEnergytime(iTime,1) = thisTime(1);
%     peakEnergytime(iTime,2) = thisTime(end);
% end

%%
% figure;
% plot(datetime(datevec(peakEnergytime(1:115,1))),peakEnergy(1:115));
% ylim([0 150]);
% xlim([datetime(datevec('26 Mar 2008 11:14')),datetime(datevec('26 Mar 2008 11:24'))]);
% title('Electron energy that is most correlated with optical emissions');
% ylabel('[keV]');
% %% Find correlation
% for i=1:1:length(energySlice)
%     r1(:,i)=xcorr(pixels2,energyFluxLog2(i,:));
% end
% [~,indx] = max(r1');
% peakEnergy = energySlice(indx);
%% Plot
figure;
hold on;
plot(energySlice,r);
ylabel('Correlation coefficient r_0_0');
xlabel('Electron Energy [keV]');
xlim([0,150]);
ylim([0,1]);
% timeCorr = time(timeMinIndx:timeMaxIndx);
% figure; plot(datetime(datevec(timeCorr)),peakEnergy(91:-1:1));
% xlim([datetime('2008-03-26 11:00'), datetime('2008-03-26 11:20')]);



function y1 = normalize_matrix(y)
    y1 = y-mean(y,2);
    y1 = y1./max(abs(y1'))';
end

