clear all;
opengl('save', 'software');
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-energyFlux.h5';
% h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20101018.001_bc_2min-energyFlux.h5';
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Temp\20080326.001_bc_15sec-energyFlux_v85.h5';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';

%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');

energySlice = 100; %keV
%%
dascImage = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
[dascKeo,dascLat,meridian] = create_keogram(dascImage,dascData.latitude,dascData.longitude);
pfisrImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
%%
pfisrImageNumF = pfisrImage./permute(repmat(pfisrData.zEnergyBin,1,1,size(pfisrImage,1)),[3,1,2]);
%%
[pfisrImage150,pfisrLat150,pfisrLon150,projAlt150] = get_pfisr_energy_slice(pfisrImage,pfisrData,150,100); 
[pfisrKeo150,pfisrPar150] = create_keogram(pfisrImage150,pfisrLat150,pfisrLon150,'meridian',meridian,'latPixelNum',6);
[pfisrImage100,pfisrLat100,pfisrLon100,projAlt100] = get_pfisr_energy_slice(pfisrImage,pfisrData,100,100); 
[pfisrKeo100,pfisrPar100] = create_keogram(pfisrImage100,pfisrLat100,pfisrLon100,'meridian',meridian,'latPixelNum',6);
[pfisrImage70,pfisrLat70,pfisrLon70,projAlt70] = get_pfisr_energy_slice(pfisrImage,pfisrData,70,100); 
[pfisrKeo70,pfisrPar70] = create_keogram(pfisrImage70,pfisrLat70,pfisrLon70,'meridian',meridian,'latPixelNum',6);
[pfisrImage30,pfisrLat30,pfisrLon30,projAlt30] = get_pfisr_energy_slice(pfisrImage,pfisrData,30,100); 
[pfisrKeo30,pfisrPar30] = create_keogram(pfisrImage30,pfisrLat30,pfisrLon30,'meridian',meridian,'latPixelNum',6);
%%
[pfisrImageAll,pfisrLatAll,pfisrLonAll,projAltAll]=integrate_pfisr_energy(pfisrImage,pfisrData,100);
[pfisrKeoAll,pfisrParAll] = create_keogram(pfisrImageAll,pfisrLatAll,pfisrLonAll,'meridian',meridian,'latPixelNum',6);
%%
[pfisrImageLE,pfisrLatLE,pfisrLonLE,projAltLE]=integrate_pfisr_energy(pfisrImage,pfisrData,15,[1,50].*1000);
[pfisrKeoLE,pfisrParLE] = create_keogram(pfisrImageLE,pfisrLatLE,pfisrLonLE,'meridian',meridian,'latPixelNum',6);

[pfisrImageHE,pfisrLatHE,pfisrLonHE,projAltHE]=integrate_pfisr_energy(pfisrImage,pfisrData,100,[50,300].*1000);
[pfisrKeoHE,pfisrParHE] = create_keogram(pfisrImageHE,pfisrLatHE,pfisrLonHE,'meridian',meridian,'latPixelNum',6);

%% pfisr diff number density
[pfisrImageLENumF,pfisrLatLENumF,pfisrLonLENumF,projAltLENumF]=integrate_pfisr_energy(pfisrImageNumF,pfisrData,15,[1,50].*1000);
[pfisrKeoLENumF,pfisrParLENumF] = create_keogram(pfisrImageLENumF,pfisrLatLENumF,pfisrLonLENumF,'meridian',meridian,'latPixelNum',6);

[pfisrImageHENumF,pfisrLatHENumF,pfisrLonHENumF,projAltHENumF]=integrate_pfisr_energy(pfisrImageNumF,pfisrData,100,[50,300].*1000);
[pfisrKeoHENumF,pfisrParHENumF] = create_keogram(pfisrImageHENumF,pfisrLatHENumF,pfisrLonHENumF,'meridian',meridian,'latPixelNum',6);

multiWaitbar('CLOSEALL');
%%
timeMinStr = '26-Mar-2008 10:30';
timeMaxStr = '26-Mar-2008 11:45';

totalPanelNo=7;
clf;

p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);

q=p(1);

q(1).select();
colormap(gca,'viridis');
ax=plot_2D_time_series(dascData.time,dascLat,dascKeo,0.25,0,timeMinStr,timeMaxStr);
% label_time_axis(time, true, 1/6,timeMinStr,timeMaxStr);
hold on; 
timeMinIndx = find_time(pfisrData.time,'26-Mar-2008 11:05');
timeMaxIndx = find_time(pfisrData.time,'26-Mar-2008 11:30');
peakLatitude150 = find_peak_latitude(pfisrKeo150,pfisrPar150);
peakLatitude100 = find_peak_latitude(pfisrKeo100,pfisrPar100);
peakLatitude70 = find_peak_latitude(pfisrKeo70,pfisrPar70);
peakLatitude30 = find_peak_latitude(pfisrKeo30,pfisrPar30);
% peakLatitude10 = find_peak_latitude(pfisrKeo10,pfisrPar10);
h150 = plot(pfisrData.time(timeMinIndx:timeMaxIndx-50),peakLatitude150(timeMinIndx:timeMaxIndx-50),'b'); %150
h100 = plot(pfisrData.time(timeMinIndx:timeMaxIndx-35),peakLatitude100(timeMinIndx:timeMaxIndx-35),'g'); %100
h70  = plot(pfisrData.time(timeMinIndx:timeMaxIndx-30),peakLatitude70(timeMinIndx:timeMaxIndx-30),'y'); %70
h30  = plot(pfisrData.time(timeMinIndx+35:timeMaxIndx-20),peakLatitude30(timeMinIndx+35:timeMaxIndx-20),'red'); %30
% plot(pfisrData.time(timeMinIndx:timeMaxIndx),peakLatitude10(timeMinIndx:timeMaxIndx),'k');
legend([h150 h100 h70 h30],'150 keV', '100 keV', '70 keV', '30 keV','Location','NorthWest');

caxis([300 400]);
axPos = get(gca, 'position');
c = colorbar('eastoutside');

set(gca,'position',axPos);

% reducing color bar thickness
cPos=get(c,'Position');
cPos(3)=0.2*cPos(3);
set(c, 'Position',cPos);
ylim([64.9 65.6]);
C = define_universal_constants();
q(2).select();
% plot_pfisr_energy_keogram(pfisrKeo150,pfisrPar150,pfisrData,false,timeMinStr,timeMaxStr,'[eV/m^2 sr s eV]');
% ylabel('150 keV');
plot_pfisr_energy_keogram(pfisrKeoAll*C.e*pi*1000,pfisrParAll,pfisrData,false,timeMinStr,timeMaxStr,'[mW/m^2]',[-2 1]);
ylabel('1-300 keV');

q(3).select();
plot(pfisrData.time,pfisrKeoLE(3,:)*C.e*pi*(1e+7).*(1e-4),'k'); 
hold on; 
plot(pfisrData.time,pfisrKeoHE(3,:)*C.e*pi*(1e+7).*(1e-4),'r');
set(gca,'XTickLabel','');
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
legend('1-50 keV','50-300 keV');
ylabel('[mW/m^2]');
ylim([0,3]);

q(4).select();
% plot(pfisrData.time,pfisrKeoLENumF(3,:)*pi.*(1e-4),'k'); 
% hold on; 
plot(pfisrData.time,pfisrKeoHENumF(3,:)*pi.*(1e-4),'r'); % convert Sr, and m2 to cm2
set(gca,'XTickLabel','');
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
legend('50-300 keV');
ylabel('F_a_t_m [cm^-^2 s^-^1]');
ylim([10^5,10^7]);

% q(5).select();
% % plot(pfisrData.time,pfisrKeoLENumF(3,:)*pi.*(1e-4),'k'); 
% % hold on; 
% plot(pfisrData.time,(pfisrKeo150(3,:)*pi.*(1e-4)./(150e+3))*70*1000,'r');
% set(gca,'XTickLabel','');
% xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
% legend('150 keV');
% ylabel('F_a_t_m [cm^-^2 s^-^1]');
% % ylim([10^5,10^7]);

q(5).select();
plot_pfisr_energy_keogram(pfisrKeo100,pfisrPar100,pfisrData,false,timeMinStr,timeMaxStr,'[eV/m^2 sr s eV]');
ylabel('100 keV');

q(6).select();
plot_pfisr_energy_keogram(pfisrKeo70,pfisrPar70,pfisrData,false,timeMinStr,timeMaxStr,'[eV/m^2 sr s eV]');
ylabel('50 keV');

q(7).select();
plot_pfisr_energy_keogram(pfisrKeo30,pfisrPar30,pfisrData,false,timeMinStr,timeMaxStr,'[eV/m^2 sr s eV]');
label_time_axis(pfisrData.time,true,1/6,timeMinStr,timeMaxStr);
ylabel('30 keV');


%%
function [image,lat,lon,projectionAltitude]=get_pfisr_energy_slice(pfisrImage,pfisrData,energySlice,projectionEnergySlice)
    for i=1:1:length(pfisrData.time)
        for ib = 1:1:size(pfisrData.latitude,1)
        image(i,ib) = interp1(pfisrData.zEnergyBin(ib,:),reshape(pfisrImage(i,ib,:),1,[]),energySlice*1000);
        end
    end
    projectionAltitude = calculate_peak_altitude_of_ionization(projectionEnergySlice*1000,pfisrData.timeNeutralAtmosphere);
    for ib = 1:1:size(pfisrData.latitude,1)
        lat(ib) = interp1(pfisrData.zEnergyBin(ib,:),pfisrData.latitude(ib,:),projectionEnergySlice*1000);
        lon(ib) = interp1(pfisrData.zEnergyBin(ib,:),pfisrData.longitude(ib,:),projectionEnergySlice*1000);
    end
end
function [image,lat,lon,projectionAltitude]=integrate_pfisr_energy(pfisrImage,...
    pfisrData,projectionEnergySlice,energyLim)
    if nargin <4
        energyLim = [1,1000].*1000;
    end
    energyMinIndx =  interp1(pfisrData.zEnergyBin(1,:),...
        1:1:length(pfisrData.zEnergyBin(1,:)),energyLim(1),'nearest');
    energyMaxIndx =  interp1(pfisrData.zEnergyBin(1,:),...
        1:1:length(pfisrData.zEnergyBin(1,:)),energyLim(2),'nearest');
    energyIndx = energyMinIndx:1:energyMaxIndx;
    for i=1:1:length(pfisrData.time)
        
    image(i,:) = trapz(pfisrData.zEnergyBin(1,energyIndx),squeeze(pfisrImage(i,:,energyIndx)),2);     
        
    end
    projectionAltitude = calculate_peak_altitude_of_ionization(projectionEnergySlice*1000,pfisrData.timeNeutralAtmosphere);
    for ib = 1:1:size(pfisrData.latitude,1)
        lat(ib) = interp1(pfisrData.zEnergyBin(ib,energyIndx),pfisrData.latitude(ib,energyIndx),projectionEnergySlice*1000);
        lon(ib) = interp1(pfisrData.zEnergyBin(ib,energyIndx),pfisrData.longitude(ib,energyIndx),projectionEnergySlice*1000);
    end
end

function plot_pfisr_energy_keogram(pfisrKeo,pfisrLat,pfisrData,setLabelTime,...
    timeMinStr,timeMaxStr,cLabelStr,cLim)
    pfisrKeo(isnan(pfisrKeo))=10^-3;
    pfisrKeo(pfisrKeo<=0) = 10^-3;
    colormap(gca,'inferno');
    plot_2D_time_series(pfisrData.time,pfisrLat,log10(pfisrKeo),0.5,0,timeMinStr,timeMaxStr);
    if setLabelTime==true
        label_time_axis(pfisrData.time, true, 1/8,timeMinStr,timeMaxStr);
    end
    if nargin < 7
        cLabelStr = '';
    end
    if nargin < 8
        cLim = [7 10];
    end
    axPos = get(gca, 'position');
    c1 = colorbar('eastoutside','Limits',cLim);
    caxis(cLim);
    set(gca,'position',axPos);
    % set(ax1,'YLim',get(ax,'YLim'));
    % reducing color bar thickness
    cPos=get(c1,'Position');
    cPos(3)=0.2*cPos(3);
    set(c1, 'Position',cPos);
    ylim([64.9 65.6]);
    ylabel(c1,cLabelStr,'FontSize',8);
    
end

function peakLatitude=find_peak_latitude(keo,lat)
    [~,index] = max(keo,[],1);
    peakLatitude = movmean(lat(index),12); 
end


