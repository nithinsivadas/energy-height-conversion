clear all;
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-energyFlux.h5';
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20101018.001_bc_2min-energyFlux.h5';
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20080326.001_bc_15sec-energyFlux_v85.h5';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';

dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');
time = dascData.time;

timeMinIndx = find_time(time,'26-Mar-2008 10:30');
timeMaxIndx = find_time(time,'26-Mar-2008 11:45');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');
%%
ibeam = 13;
energySlice = pfisrData.zEnergyBin(ibeam,:)/1000; %keV

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
timeMinIndx = find_time(thisTimeDASC,'26-Mar-2008 10:30');
% timeMaxIndx = find_time(thisTimeDASC,'26-Mar-2008 11:22');
dTimeIndx = 16;
for iTime = 1:1:length(thisTimeDASC)-dTimeIndx-1
%     iTime=1;
    iTimeIndx = timeMinIndx+iTime-1:1:timeMinIndx+dTimeIndx+iTime-1;
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
    rArray(iTime,:) = r;
    [value,indx] = max(r);
    if value<0.01
        peakEnergy(iTime) = nan;
    else
        peakEnergy(iTime) = energySlice(indx);
    end
    peakEnergytime(iTime,1) = thisTime(1);
    peakEnergytime(iTime,2) = thisTime(end);
end

%%
totalPanelNo=2;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);
q=p(1);
q(1).select();

% 4e. Correlation 
timeMinStr = '26-Mar-2008 10:30:00';
timeMaxStr = '26-Mar-2008 11:40:00';
tindx1 = find_time(peakEnergytime(:,1),timeMinStr);
tindx2 = find_time(peakEnergytime(:,1),timeMaxStr);
sz=max(rArray(tindx1:tindx2,:)');
sz(sz<=0.01)=0.01;

yyaxis right    
area(peakEnergytime(tindx1:tindx2,1),...
    sz,'FaceColor',[0.8,0.8,0.8],'EdgeColor',[0.8,0.8,0.8],...
    'FaceAlpha',0.5,'EdgeAlpha',0.5);
ylim([0 1]);
ylabel({'Correlation','coefficient','[a.u.]'});
% grid on;
set(gca,'ycolor',[0.5,0.5,0.5],'YTick',[0,0.5,1],'Ygrid','on');
label_time_axis(pfisrData.time,false,1/6,timeMinStr,timeMaxStr);
hold on;
line([peakEnergytime(tindx1,1),peakEnergytime(tindx2)],...
    [0.5,0.5],'Color',[0.8,0.8,0.8]);
yyaxis left
scatter(peakEnergytime(tindx1:tindx2,1)...
    ,peakEnergy(tindx1:tindx2),(sz*5).^3,...
    'filled','MarkerFaceColor','r', 'MarkerFaceAlpha',.5);
set(gca,'ycolor','r','YTick',[0,30,50,100,150],'Ygrid','on')
label_time_axis(pfisrData.time,false,1/6,timeMinStr,timeMaxStr);
ylim([0 150]);
ylabel({'Max e^- energy','that correlates','with emissions','[keV]'});

% 4e correlation with
q(2).select();
plot_2D_time_series(peakEnergytime(:,1),energySlice,real(rArray)',0.5,0,timeMinStr,timeMaxStr);
label_time_axis(pfisrData.time,true,1/6,timeMinStr,timeMaxStr);
% colormap(viridis);
colormap(get_colormap('w','r'));
set(gca,'YScale','log','YTick',[1,10,30,100,300,1000],'YGrid','on','YMinorGrid','off');
caxis([0 +1]);
ylabel({'e^- energy flux','correlation','with emissions','[keV]'});

% reducing color bar thickness
c=colorbar_thin('Width',0.2,'YLabel','r_0_0 [a.u.]');
% axPos = get(gca, 'position');
% c = colorbar('eastoutside');
% set(gca,'position',axPos);
% 
% cPos=get(c,'Position');
% cPos(3)=0.2*cPos(3);
% set(c, 'Position',cPos);
% ylabel(c,'r_0_0 [a.u.]')
scatter(peakEnergytime(tindx1:tindx2,1)...
    ,peakEnergy(tindx1:tindx2),8,...
    'filled','MarkerFaceColor','k', 'MarkerFaceAlpha',.4);

% %% Find correlation
% for i=1:1:length(energySlice)
%     r1(:,i)=xcorr(pixels2,energyFluxLog2(i,:));
% end
% [~,indx] = max(r1');
% peakEnergy = energySlice(indx);
%% Plot
figure;
hold on;
indx1 = find_time(peakEnergytime(:,1),'26-Mar-2008 11:16:18');
plot(energySlice,rArray(indx1,:));
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

