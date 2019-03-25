%% Paper 2: Figure 4 ver 2
clear all;
opengl('save', 'software');
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
h5FileStrDASC = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
%%
mspFileStr = 'C:\Users\nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 2\Data\msp_vvel.mat';
dataMSP = get_msp_data(mspFileStr,'26 Mar 2008 10:00','26 Mar 2008 12:00');
%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStrDASC,...
    'plotModeStr','OpticalImage');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');

energySlice = 100; %keV

%% 
dascImage = permute(h5read(h5FileStrDASC,'/DASC/ASI'),[3 2 1]);
%% Cropping DASC
minTimeStr = '26 Mar 2008 10:30';
maxTimeStr = '26 Mar 2008 12:00';
timeMinIndx = find_time(dascData.time,minTimeStr);
timeMaxIndx = find_time(dascData.time,maxTimeStr);
dascTime  = dascData.time(timeMinIndx:timeMaxIndx);
dascImage1 = dascImage(timeMinIndx:timeMaxIndx,:,:);

%% Keogram DASC
[dascKeo,dascLat,meridian] = create_keogram(dascImage1,dascData.latitude,dascData.longitude,...
    'time',datenum('26 Mar 2008 11:00'));
%% PFISR Keogram
pfisrImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
pfisrImageNumF = pfisrImage./permute(repmat(pfisrData.zEnergyBin,1,1,size(pfisrImage,1)),[3,1,2]);
%%
% timeMinIndx = find_time(dascData.time,minTimeStr);
% timeMaxIndx = find_time(dascData.time,maxTimeStr);
% pfisrImage1 = pfisrImage(timeMinIndx:timeMaxIndx,:,:);
% pfisrImageNumF1 = pfisrImageNumF(timeMinIndx:timeMaxIndx,:,:);
%%
% E = 200;
% projAlt = pfisrData.projectionAltitude(find_altitude(pfisrData.zEnergyBin(1,:),E*1000));
% [pfisrImage150,pfisrLat150,pfisrLon150,projAlt150] = get_pfisr_energy_slice(pfisrImage,pfisrData,E,projAlt); 
% [pfisrKeo150,pfisrPar150] = create_keogram(pfisrImage150,pfisrLat150,pfisrLon150,'latPixelNum',6);
%%
E = 100;
projAlt = pfisrData.projectionAltitude(find_altitude(pfisrData.zEnergyBin(1,:),E*1000));
[pfisrImage100,pfisrLat100,pfisrLon100,projAlt100] = get_pfisr_energy_slice(pfisrImage,pfisrData,E,projAlt); 
[pfisrKeo100,pfisrPar100] = create_keogram(pfisrImage100,pfisrLat100,pfisrLon100,'latPixelNum',6);

E = 70;
projAlt = pfisrData.projectionAltitude(find_altitude(pfisrData.zEnergyBin(1,:),E*1000));
[pfisrImage70,pfisrLat70,pfisrLon70,projAlt70] = get_pfisr_energy_slice(pfisrImage,pfisrData,E,projAlt); 
[pfisrKeo70,pfisrPar70] = create_keogram(pfisrImage70,pfisrLat70,pfisrLon70,'latPixelNum',6);

E = 30;
projAlt = pfisrData.projectionAltitude(find_altitude(pfisrData.zEnergyBin(1,:),E*1000));
[pfisrImage30,pfisrLat30,pfisrLon30,projAlt30] = get_pfisr_energy_slice(pfisrImage,pfisrData,E,projAlt); 
[pfisrKeo30,pfisrPar30] = create_keogram(pfisrImage30,pfisrLat30,pfisrLon30,'latPixelNum',6);

%%
% meridian = convert_longitude(meridian);
% [pfisrImage150,pfisrLat150,pfisrLon150,projAlt150] = get_pfisr_energy_slice(pfisrImage,pfisrData,150,100); 
% [pfisrKeo150,pfisrPar150] = create_keogram(pfisrImage150,pfisrLat150,pfisrLon150,'meridian',meridian,'latPixelNum',6);
% [pfisrImage100,pfisrLat100,pfisrLon100,projAlt100] = get_pfisr_energy_slice(pfisrImage,pfisrData,100,100); 
% [pfisrKeo100,pfisrPar100] = create_keogram(pfisrImage100,pfisrLat100,pfisrLon100,'meridian',meridian,'latPixelNum',6);
% [pfisrImage70,pfisrLat70,pfisrLon70,projAlt70] = get_pfisr_energy_slice(pfisrImage,pfisrData,70,100); 
% [pfisrKeo70,pfisrPar70] = create_keogram(pfisrImage70,pfisrLat70,pfisrLon70,'meridian',meridian,'latPixelNum',6);
% [pfisrImage30,pfisrLat30,pfisrLon30,projAlt30] = get_pfisr_energy_slice(pfisrImage,pfisrData,30,100); 
% [pfisrKeo30,pfisrPar30] = create_keogram(pfisrImage30,pfisrLat30,pfisrLon30,'meridian',meridian,'latPixelNum',6);
%%
[pfisrImageAll,pfisrLatAll,pfisrLonAll,projAltAll]=integrate_pfisr_energy(pfisrImage,pfisrData,110);
[pfisrKeoAll,pfisrParAll] = create_keogram(pfisrImageAll,pfisrLatAll,pfisrLonAll,'latPixelNum',6);
%%
[pfisrImageLE,pfisrLatLE,pfisrLonLE,projAltLE]=integrate_pfisr_energy(pfisrImage,pfisrData,5,[1,10].*1000);
[pfisrKeoLE,pfisrParLE] = create_keogram(pfisrImageLE,pfisrLatLE,pfisrLonLE,'latPixelNum',6);
%%
[pfisrImageME1,pfisrLatME1,pfisrLonME1,projAltME1]=integrate_pfisr_energy(pfisrImage,pfisrData,20,[10,30].*1000);
[pfisrKeoME1,pfisrParME1] = create_keogram(pfisrImageME1,pfisrLatME1,pfisrLonME1,'latPixelNum',6);

%% 
[pfisrImageME,pfisrLatME,pfisrLonME,projAltME]=integrate_pfisr_energy(pfisrImage,pfisrData,70,[30,100].*1000);
[pfisrKeoME,pfisrParME] = create_keogram(pfisrImageME,pfisrLatME,pfisrLonME,'latPixelNum',6);
%%
[pfisrImageHE,pfisrLatHE,pfisrLonHE,projAltHE]=integrate_pfisr_energy(pfisrImage,pfisrData,150,[100,300].*1000);
[pfisrKeoHE,pfisrParHE] = create_keogram(pfisrImageHE,pfisrLatHE,pfisrLonHE,'latPixelNum',6);

%% pfisr diff number density
% [pfisrImageLENumF,pfisrLatLENumF,pfisrLonLENumF,projAltLENumF]=integrate_pfisr_energy(pfisrImageNumF,pfisrData,15,[5,30].*1000);
% [pfisrKeoLENumF,pfisrParLENumF] = create_keogram(pfisrImageLENumF,pfisrLatLENumF,pfisrLonLENumF,'latPixelNum',6);
% 
% [pfisrImageMENumF,pfisrLatMENumF,pfisrLonMENumF,projAltMENumF]=integrate_pfisr_energy(pfisrImageNumF,pfisrData,70,[50,1000].*1000);
% [pfisrKeoMENumF,pfisrParMENumF] = create_keogram(pfisrImageMENumF,pfisrLatMENumF,pfisrLonMENumF,'latPixelNum',6);
% 
[pfisrImageHENumF,pfisrLatHENumF,pfisrLonHENumF,projAltHENumF]=integrate_pfisr_energy(pfisrImageNumF,pfisrData,150,[100,300].*1000);
[pfisrKeoHENumF,pfisrParHENumF] = create_keogram(pfisrImageHENumF,pfisrLatHENumF,pfisrLonHENumF,'latPixelNum',6);

multiWaitbar('CLOSEALL');
%%
timeMinStr = '26-Mar-2008 10:30';
timeMaxStr = '26-Mar-2008 11:40';

totalPanelNo=3;
% clf;

p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);

q=p(1);

% 4a. Keogram, White light image
q(1).select();
% colormap(gca,get_colormap('k',[0,1,0]));
colormap(gca,viridis);

ax=plot_2D_time_series(dascTime,dascLat,dascKeo,0.25,0,timeMinStr,timeMaxStr);
hold on; 
timeMinIndx = find_time(pfisrData.time,'26-Mar-2008 11:15');
timeMaxIndx = find_time(pfisrData.time,'26-Mar-2008 11:24');
% peakLatitude150 = find_peak_latitude(pfisrKeo150,pfisrPar150);
peakLatitude100 = find_peak_latitude(pfisrKeo100,pfisrPar100);
peakLatitude70 = find_peak_latitude(pfisrKeo70,pfisrPar70);
peakLatitude30 = find_peak_latitude(pfisrKeo30,pfisrPar30);

% h150 = plot(pfisrData.time(timeMinIndx:timeMaxIndx-50),peakLatitude150(timeMinIndx:timeMaxIndx-50),'b','LineWidth',1); %150
h100 = plot(pfisrData.time(timeMinIndx-10:timeMaxIndx-15),peakLatitude100(timeMinIndx-10:timeMaxIndx-15),'c','LineWidth',1); %100
h70  = plot(pfisrData.time(timeMinIndx-5:timeMaxIndx-10),peakLatitude70(timeMinIndx-5:timeMaxIndx-10),'Color',[0.9,0.9,1],'LineWidth',1); %70
h30  = plot(pfisrData.time(timeMinIndx+8:timeMaxIndx),peakLatitude30(timeMinIndx+8:timeMaxIndx),'Color',[1,0.7,0.7],'LineWidth',1); %30

legend([h100 h70 h30], '100 keV', '70 keV', '30 keV','Location','NorthWest');
% caxis([0.15 0.35]);
caxis([350 380]);
axPos = get(gca, 'position');
c = colorbar('eastoutside');
set(gca,'position',axPos);

% reducing color bar thickness
cPos=get(c,'Position');
cPos(3)=0.2*cPos(3);
set(c, 'Position',cPos);
ylim([64.8 65.6]);
% ylim([64.7 65.6]);

ylabel({'White light','emission',['lat [N',char(176),']']});
label_time_axis(pfisrData.time,false,1/6,timeMinStr,timeMaxStr);

% Figure 4b 100 keV Electron Energy Flux
q(2).select();
plot_pfisr_energy_keogram(pfisrKeo100,pfisrPar100,pfisrData,false,timeMinStr,timeMaxStr,'[eV/m^2 sr s eV]');
ylabel({'e^- Energy Flux','100 keV',['lat [N',char(176),']']});
label_time_axis(pfisrData.time,false,1/6,timeMinStr,timeMaxStr);

% % Figure 4c 1-300 keV Total energy flux 
% q(3).select();
C = define_universal_constants();
% plot_pfisr_energy_keogram(pfisrKeoAll*C.e*pi*1000,pfisrParAll,pfisrData,false,timeMinStr,timeMaxStr,'[mW/m^2]',[-1 1]);
% ylabel({'e^- Energy Flux','1-300 keV','lat [N^0]'});
% label_time_axis(pfisrData.time,false,1/6,timeMinStr,timeMaxStr);

q(3).select();
x = pfisrData.time;
iLat = 5;
y = [pfisrKeoHE(iLat,:);...
    pfisrKeoME(iLat,:)-pfisrKeoHE(iLat,:);...
    pfisrKeoME1(iLat,:);...
    pfisrKeoLE(iLat,:)].*C.e*pi*(1e+7).*(1e-4); % It is pi and not 2pi, cause you are only looking at downgoing electrons
dascLatIndx = find_altitude(dascLat,pfisrParHE(iLat));

b = area(x,y');
b(1).FaceColor = 'm';
b(2).FaceColor = 'y';
b(3).FaceColor = [0.5, 0.8, 0];
b(4).FaceColor = [0.5, 0.5, 0.5];
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(3).EdgeColor = 'none';
b(4).EdgeColor = 'none';
set(gca,'XTickLabel','');

% yyaxis right;
% a=plot(dascTime, dascKeo(dascLatIndx,:),'b');
% ylim([350,380]);

% yyaxis left;
xlim([datenum(timeMinStr), datenum(timeMaxStr)]);
% legend([b(4),b(3),b(2),b(1),a],{'1-10','10-30 keV','30-50 keV','50-300 keV','Intensity'},'Location','NorthWest');
legend([b(4),b(3),b(2),b(1)],{'1-10','10-30 keV','30-100 keV','100-300 keV'},'Location','NorthWest');
ylabel({'e^- Energy Flux','[mW/m^2]'});
ylim([0,2]);
label_time_axis(pfisrData.time,true,1/6,timeMinStr,timeMaxStr);

%% Plot the net power of electron > 30 keV
figure; 
plot(pfisrData.time,(pfisrKeoHE(5,:)+pfisrKeoME(5,:)).*C.e*pi*(1e+7).*(1e-4)); 
xlim([datenum(timeMinStr), datenum(timeMaxStr)]); 
label_time_axis(pfisrData.time,true,1/6,timeMinStr,timeMaxStr); 
ylim([0 1]);
title('Energy Flux >30 keV');
ylabel('Power [mW/m^2]');

%% Plot the net diff flux of electron > 100 keV
figure; 
plot(pfisrData.time,(pfisrKeoHENumF(5,:)).*pi*(1e-4)); 
xlim([datenum(timeMinStr), datenum(timeMaxStr)]); 
label_time_axis(pfisrData.time,true,1/6,timeMinStr,timeMaxStr); 
% ylim([0 1]);
title('Differential Number Flux >100 keV');
ylabel('[m-2 s-1]');


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

function [glat,glon,galt] = convert_azel_to_glatlon(az,el,galt,glat0,glon0,h0,time0)
    C = define_universal_constants;
%     slantRange = -C.RE.*sind(el) + modSign(el).*sqrt((C.RE^2).*(sind(el)).^2 + galt.*1000.*(galt.*1000+2.*C.RE)); 
    slantRange = galt.*1000./sind(el);
    slantRange(isinf(slantRange)) =slantRange(2);
% slant Range = RE*sin(el) +/- sqrt(RE^2 sin^2(el) + H(H+2RE); RE - radius
% of earth, H- altitude, el - elevation
    [glat,glon,galt] = aer2geodetic(az,el,slantRange./1000,glat0,glon0,h0,wgs84Ellipsoid('km'));    
end


