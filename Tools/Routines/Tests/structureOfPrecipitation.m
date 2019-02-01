%% Figure 1e (Variation of structure)
clear all;
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
h5FileStrNe = 'G:\My Drive\Research\Projects\Paper 2\Data\Event 1\20080326.001_bc_15sec-fitcal.h5';
%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');

dascImage = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
pfisrEnergyFluxImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
pfisrNeImage = permute(h5read(h5FileStr,'/inputData/Ne'),[3 2 1]);

%% Trim data
timeMinStr = '26 Mar 2008 10:30';
timeMaxStr = '26 Mar 2008 11:40';
%%
asi.indx = (find_time(dascData.time,timeMinStr):1:find_time(dascData.time,timeMaxStr))';
asi.data = dascImage(asi.indx,:,:);
asi.lat = permute(reshape(repmat(dascData.latitude(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.lon = permute(reshape(repmat(dascData.longitude(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.az = permute(reshape(repmat(dascData.azimuth(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.el = permute(reshape(repmat(dascData.elevation(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.time = repmat(dascData.time(asi.indx),1,size(asi.data,2),size(asi.data,3));
%%
pfisr.indx = (find_time(pfisrData.time,timeMinStr):1:find_time(pfisrData.time,timeMaxStr))';
pfisr.energyFlux = pfisrEnergyFluxImage(pfisr.indx,:,:);
pfisr.lat = permute(reshape(repmat(pfisrData.latitude(:),1,length(pfisr.indx)),size(pfisr.energyFlux,2),size(pfisr.energyFlux,3),length(pfisr.indx)),[3,1,2]);
pfisr.lon = permute(reshape(repmat(pfisrData.longitude(:),1,length(pfisr.indx)),size(pfisr.energyFlux,2),size(pfisr.energyFlux,3),length(pfisr.indx)),[3,1,2]);
pfisr.zEnergyBin = permute(reshape(repmat(pfisrData.zEnergyBin(:),1,length(pfisr.indx)),size(pfisr.energyFlux,2),size(pfisr.energyFlux,3),length(pfisr.indx)),[3,1,2]);
pfisr.time = repmat(pfisrData.time(pfisr.indx),1,size(pfisr.energyFlux,2),size(pfisr.energyFlux,3));

pfisrNe.indx = pfisr.indx;
pfisrNe.data = pfisrNeImage(pfisrNe.indx,:,:);
tempCoords = permute(h5read(h5FileStr,'/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt'),[3 2 1]);
pfisrNe.lat =  repmat(tempCoords(1,:,:),length(pfisrNe.indx),1,1);
pfisrNe.lon = repmat(tempCoords(2,:,:),length(pfisrNe.indx),1,1);
pfisrNe.alt = repmat(tempCoords(3,:,:),length(pfisrNe.indx),1,1);

%%
amisr = read_amisr(h5FileStrNe);
%%
amisrData = aer_to_field_aligned_coords(amisr,60);
%%
% temp.electronDensity = amisrData.electronDensity;
% temp.xEast = repmat(amisrData.cartCoords.xEast,1,1,length(amisrData.time(1,:)));
% temp.yNorth = repmat(amisrData.cartCoords.yNorth,1,1,length(amisrData.time(1,:)));
% temp.zUp = repmat(amisrData.cartCoords.zUp,1,1,length(amisrData.time(1,:)));
% temp.xEast0 = repmat(amisrData.origCartCoords.xEast,1,1,length(amisrData.time(1,:)));
% temp.yNorth0 = repmat(amisrData.origCartCoords.yNorth,1,1,length(amisrData.time(1,:)));
% temp.zUp0 = repmat(amisrData.origCartCoords.zUp,1,1,length(amisrData.time(1,:)));

%%
% pfisrNe0.indx = (find_time(amisrData.time(1,:),timeMinStr):1:find_time(amisrData.time(1,:),timeMaxStr))';
% pfisrNe0.time = amisrData.time(1,pfisrNe0.indx)';
% multiWaitbar('Interpolating...',0);
% id = 1./length(pfisrNe0.time);
% for iT = 1:1:length(pfisrNe0.time)
%     ne = squeeze(amisrData.electronDensity(:,:,pfisrNe0.indx(iT)));
%     F3 = scatteredInterpolant(amisrData.origCartCoords.xEast(:),...
%                               amisrData.origCartCoords.yNorth(:),...
%                               amisrData.origCartCoords.zUp(:),...
%                               ne(:),'natural');
%     pfisrNe0.data(:,:,iT) = F3(amisrData.cartCoords.xEast,...
%                                 amisrData.cartCoords.yNorth,...
%                                 amisrData.cartCoords.zUp); 
%     multiWaitbar('Interpolating...','Increment',id);
% end
% 
% [pfisrNe0.lat,pfisrNe0.lon,pfisrNe0.alt]=...
%     ned2geodetic(amisrData.cartCoords.yNorth,amisrData.cartCoords.xEast,...
%     -amisrData.cartCoords.zUp,amisrData.site.latitude,amisrData.site.longitude,...
%     amisrData.site.altitude./1000,wgs84Ellipsoid('km'));
% 
% pfisrNe0.data = permute(pfisrNe0.data,[3,2,1]);
% pfisrNe0.lat = pfisrNe0.lat';
% pfisrNe0.lon = pfisrNe0.lon';
% pfisrNe0.alt = pfisrNe0.alt';

%%
% pfisrNe0.indx = (find_time(amisrData.time(1,:),timeMinStr):1:find_time(amisrData.time(1,:),timeMaxStr))';
% pfisrNe0.time = amisrData.time(1,pfisrNe0.indx)';
% multiWaitbar('Interpolating...',0);
% id = 1./length(pfisrNe0.time);
% for iT = 1:1:length(pfisrNe0.time)
%     ne = squeeze(amisrData.electronDensity(:,:,pfisrNe0.indx(iT)));
%     F3 = scatteredInterpolant(amisrData.origCartCoords.xEast(:),...
%                               amisrData.origCartCoords.yNorth(:),...
%                               amisrData.origCartCoords.zUp(:),...
%                               ne(:),'natural');
%     pfisrNe0.data(:,:,iT) = F3(amisrData.cartCoords.xEast,...
%                                 amisrData.cartCoords.yNorth,...
%                                 amisrData.cartCoords.zUp); 
%     multiWaitbar('Interpolating...','Increment',id);
% end
% 
% [pfisrNe0.lat,pfisrNe0.lon,pfisrNe0.alt]=...
%     ned2geodetic(amisrData.cartCoords.yNorth,amisrData.cartCoords.xEast,...
%     -amisrData.cartCoords.zUp,amisrData.site.latitude,amisrData.site.longitude,...
%     amisrData.site.altitude./1000,wgs84Ellipsoid('km'));
% 
% pfisrNe0.data = permute(pfisrNe0.data,[3,2,1]);
% pfisrNe0.lat = pfisrNe0.lat';
% pfisrNe0.lon = pfisrNe0.lon';
% pfisrNe0.alt = pfisrNe0.alt';
%%
% altMin =60;
% altMax = 200;
% altIndx = (find_altitude(pfisrNe0.alt(13,:),altMin):1:find_altitude(pfisrNe0.alt(13,:),altMax))';
% pfisrNe0.data = pfisrNe0.data(:,:,altIndx);
% pfisrNe0.lat = pfisrNe0.lat(:,altIndx);
% pfisrNe0.lon = pfisrNe0.lon(:,altIndx);
% pfisrNe0.alt = pfisrNe0.alt(:,altIndx);
%% Undersampled asi measurements
[I, J] = ndgrid(1:length(asi.lat),1:length(asi.lat));
asilat = squeeze(asi.lat(1,:,:));
asilon = squeeze(asi.lon(1,:,:));
asiel = squeeze(asi.el(1,:,:));
asiaz = squeeze(asi.az(1,:,:));
elIndx = asiel>10;
Fi = scatteredInterpolant(asilat(elIndx),asilon(elIndx),I(elIndx),'nearest');
Fj = scatteredInterpolant(asilat(elIndx),asilon(elIndx),J(elIndx),'nearest');
%% Pixels corresponding to radar beam positions for different energy slice & altitude slice
vi_eflux = squeeze(Fi(pfisr.lat(1,:,:),pfisr.lon(1,:,:)));
vj_eflux = squeeze(Fj(pfisr.lat(1,:,:),pfisr.lon(1,:,:)));
vi_ne = squeeze(Fi(pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:)));
vj_ne = squeeze(Fj(pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:)));
%%
% vi_ne0 = squeeze(Fi(pfisrNe0.lat,pfisrNe0.lon));
% vj_ne0 = squeeze(Fj(pfisrNe0.lat,pfisrNe0.lon));
%% Pixels within radar beam; beamWidth = 1
[i_eflux,j_eflux] = get_neighbouring_pixels(vi_eflux,vj_eflux,asiaz,asiel,1);
[i_ne,j_ne] = get_neighbouring_pixels(vi_ne,vj_ne,asiaz,asiel,1);
%%
% [i_ne0,j_ne0] = get_neighbouring_pixels(vi_ne0,vj_ne0,asiaz,asiel,1);
%% Get undersampled asi images
asi.eflux = get_undersampled_asi(asi.data,i_eflux,j_eflux);
asi.ne = get_undersampled_asi(asi.data,i_ne,j_ne);
%%
% asi.ne0 = get_undersampled_asi(asi.data,i_ne0,j_ne0);
%%
% asi.ne0 = asi.ne0(:,:,altIndx);
%% Figure of Energy Flux 
h1 = figure;
resize_figure(h1,125,200);
timeArrayStr = string(['26 Mar 2008 11:04:00';...
                       '26 Mar 2008 11:08:00';...
                       '26 Mar 2008 11:16:00';...
                       '26 Mar 2008 11:24:00';...
                       '26 Mar 2008 11:28:00']);
zValueBinArray = [3,8,30,100];
cLimArray = [9,11;...
            9,11;...
            9,11;...
            8,10];
zUnitStr = 'keV';
colorUnitStr = '[eV/m^2 sr s eV]';
titleStr = 'Energy Flux from PFISR';
generate_panels(h1,log10(pfisr.energyFlux),pfisr.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
    pfisr.zEnergyBin(1,:,:)/1000,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
colormap('inferno');

%% Electron density (descrepancy between the coordinates of electron density and energy flux)
h2 = figure;
resize_figure(h2,125,200);
timeArrayStr = string(['26 Mar 2008 11:04:00';...
                       '26 Mar 2008 11:08:00';...
                       '26 Mar 2008 11:16:00';...
                       '26 Mar 2008 11:24:00';...
                       '26 Mar 2008 11:28:00']);
zValueBinArray = [120,105,95,85];
cLimArray = [10,12;...
            10,12;...
            10,12;...
            10,12];
zUnitStr = 'km';
colorUnitStr = '[m^-^3]';
% uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1).*1000;
titleStr = 'Ne from PFISR';
pfisrNe.data(pfisrNe.data<=0) = nan;
generate_panels(h2,log10(pfisrNe.data),pfisr.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
    pfisrNe.alt(1,:,:),timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
colormap(magma);
%% Electron density (descrepancy between the coordinates of electron density and energy flux)
% h2 = figure;
% resize_figure(h2,125,200);
% timeArrayStr = string(['26 Mar 2008 11:04:00';...
%                        '26 Mar 2008 11:08:00';...
%                        '26 Mar 2008 11:16:00';...
%                        '26 Mar 2008 11:24:00';...
%                        '26 Mar 2008 11:28:00']);
% zValueBinArray = [120,105,95,85];
% cLimArray = [10,12;...
%             10,12;...
%             10,12;...
%             10,12];
% zUnitStr = 'km';
% colorUnitStr = '[m^-^3]';
% titleStr = 'Ne from PFISR';
% uniformAlt = repmat(mean(pfisrNe0.alt,1),size(pfisrNe0.alt,1),1);
% generate_panels(h2,pfisrNe0.data,pfisrNe0.time,pfisrNe0.lat,pfisrNe0.lon,...
%     uniformAlt*1000,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,...
%     titleStr,[64.6,66.2],[-149,-146]);
% colormap(magma);

%%
figure;
plot_slice(pfisrNe0.data,pfisrNe0.time,pfisrNe0.lat,pfisrNe0.lon,...
    pfisrNe0.alt*1000,'26 Mar 2008 11:00',100);
colormap(inferno);
colorbar; 

%% Figure of asi undersampled/ energy flux
h3 = figure;
resize_figure(h3,125,200);
timeArrayStr = string(['26 Mar 2008 11:04:00';...
                       '26 Mar 2008 11:08:00';...
                       '26 Mar 2008 11:16:00';...
                       '26 Mar 2008 11:24:00';...
                       '26 Mar 2008 11:28:00']);
zValueBinArray = [3,8,30,100];
cLimArray = [350,370;...
            350,370;...
            350,370;...
            350,370];
zUnitStr = 'keV';
colorUnitStr = '[a.u.]';
titleStr = 'ASI undersampled at PFISR energy flux coords';
generate_panels(h3,asi.eflux,asi.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
    pfisrData.zEnergyBin,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
colormap(viridis);

%%
figure; scatter(pfisr.lat(1,:,1),pfisr.lon(1,:,1));
hold on; 
scatter(pfisrNe.lat(1,:,120),pfisrNe.lon(1,:,120));
xlim([64.8,65.5]);
ylim([-148.2 -146.8]);
%% Figure of asi undersampled/ Ne
h4 = figure;
resize_figure(h4,125,200);
timeArrayStr = string(['26 Mar 2008 11:04:00';...
                       '26 Mar 2008 11:08:00';...
                       '26 Mar 2008 11:16:00';...
                       '26 Mar 2008 11:24:00';...
                       '26 Mar 2008 11:28:00']);
zValueBinArray = [120,105,95,85];
cLimArray = [350,370;...
            350,370;...
            350,370;...
            350,370];
zUnitStr = 'km';
colorUnitStr = '[a.u.]';
titleStr = 'ASI undersampled at PFISR Ne coords';
uniformAlt = repmat((pfisrNe.alt(1,13,:)),1,size(pfisrNe.data,2),1).*1000;
generate_panels_asi(h4,asi.ne,asi.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
    uniformAlt,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);

colormap(viridis);

%% Testing
figure;
% elIndx = asi.el(1,:,:)>30;
i = 170;
plot_DASC_aer_v1(squeeze(asi.data(i,:,:)),squeeze(asi.az(1,:,:)),squeeze(asi.el(1,:,:)));
hold on;
plot_scatter_aer(diag(asiaz(vi_eflux(:,1),vj_eflux(:,1))),diag(asiel(vi_eflux(:,1),vj_eflux(:,1))));
hold on;
plot_scatter_aer(diag(asiaz(vi_eflux(:,30),vj_eflux(:,30))),diag(asiel(vi_eflux(:,30),vj_eflux(:,30))),'c');
hold on;
linearInd = calculate_linear_index (size(asiaz), i_eflux, j_eflux, 10, 1);
% iE = 1;
% % linearInd = [];
% for iB = 1:26
%     linearInd = [linearInd;squeeze(sub2ind(size(asiaz),i_eflux(iB,iE,i_eflux(iB,iE,:)>0),j_eflux(iB,iE,j_eflux(iB,iE,:)>0)))];
% end
plot_scatter_aer(asiaz(linearInd),asiel(linearInd),'y');
plot_grid_aer(90,30,'r');
% iE = 1;
% plot_scatter_aer(diag(asiaz(i_eflux(iB,iE,i_eflux(iB,iE,:)>0),j_eflux(iB,iE,j_eflux(iB,iE,:)>0))),diag(asiel(i_eflux(iB,iE,i_eflux(iB,iE,:)>0),j_eflux(iB,iE,j_eflux(iB,iE,:)>0))),'y');
% scatter(90*vi_eflux(i,:)./1024,90*vj_eflux(i,:)./1024);
title(datestr(asi.time(i)));
%%
% figure;
% plot_DASC_aer_v1(squeeze(asi.az(1,:,:)),squeeze(asi.az(1,:,:)),squeeze(asi.el(1,:,:)));

figure;
% plot_2D_energy_slice_geodetic_v2018(pfisr.energyFlux(1,:,:),...
%     pfisr.lat(1,:,:),pfisr.lon(1,:,:),pfisrData.zEnergyBin,...
%     [],10,2,4,false,false,256);
plot_slice(pfisr.energyFlux,pfisr.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
    pfisrData.zEnergyBin,'26 Mar 2008 11:00',5);
colormap(inferno);
colorbar; 

%% Find Correlation 
[r_eflux,time_r_eflux, time2_r_eflux] = calculate_correlation_v2(pfisr.energyFlux,pfisr.time(:,1,1),asi.eflux,asi.time(:,1,1));
[r_ne,time_r_ne, time2_r_ne] = calculate_correlation_v2(pfisrNe.data,pfisr.time(:,1,1),asi.ne,asi.time(:,1,1));
%%
altIndx = find_altitude(pfisrNe.alt(1,13,:),80);
pfisrne80km = repmat(interp_nans(pfisrNe.data(:,:,altIndx)')',1,1,size(pfisrNe.data,3));
[r_ne_ne,time_r_ne_ne] = calculate_correlation_v2(pfisrNe.data,pfisr.time(:,1,1),pfisrne80km,pfisr.time(:,1,1));

%%
% [r_ne0,time_r_ne0] = calculate_correlation(pfisrNe0.data,pfisrNe0.time,asi.ne0,asi.time(:,1,1));
%%
totalPanelNo=2;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);
q=p(1);

q(1).select();
plot_2D_time_series(time_r_eflux,pfisrData.zEnergyBin(1,:)./1000,r_eflux',10/60,-1);
label_time_axis(time_r_eflux, false, 10/60);
set(gca,'YScale','log','YTick',[10,30,100,300],'YTickLabel',[10,30,100,300]);
ylabel('[keV]');
colormap(gca,get_colormap('w','r'));
grid on;
% colormap(jet);
colorbar_thin('YLabel',{'r_0','Eflux & Optical'});
caxis([0,1]);


q(2).select();
uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1);
plot_2D_time_series(time_r_ne,uniformAlt(1,1,:),r_ne',0.1,-1);
label_time_axis(time_r_ne, true, 10/60);
set(gca,'YScale','log','YTick',[70,80,100,120,200]);
colormap(gca,get_colormap('w','k'));
% title('Correlation of Ne (mag-field-aligned) with ASI');
ylabel('[Km]');
colorbar_thin('YLabel',{'r_0 (N_e & Optical)'});
caxis([0,1]);

%%
totalPanelNo=1;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);
q=p(1);

q(1).select();
uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1);
plot_2D_time_series(time_r_ne_ne,uniformAlt(1,1,:),r_ne_ne',0.1,-1);
label_time_axis(time_r_ne_ne, true, 10/60);
set(gca,'YScale','log','YTick',[70,80,100,120,200]);
colormap(gca,get_colormap('w','k'));
% title('Correlation of Ne (mag-field-aligned) with ASI');
ylabel('[Km]');
colorbar_thin('YLabel',{'r_0 (N_e@80km & N_e)'});
caxis([0,1]);
%%
figure;
plot(time_r_ne,(time2_r_ne(:)-time_r_ne(:))*24*3600);
label_time_axis(time_r_ne, true, 10/60);
% figure;
% plot_2D_time_series(time_r_ne0,mean(pfisrNe0.alt,1),r_ne0',0.1,-1);
% label_time_axis(time_r_ne0, true, 0.1);
% set(gca,'YScale','log');
% colormap(get_colormap('w','k'));
% title('Correlation of Ne with ASI');
% ylabel('[Km]');
% colorbar;
% caxis([0,1]);

%% Test plotting slice
figure;
scatter(pfisrNe.lat(1,:,86),pfisrNe.lon(1,:,86)); hold on; scatter(pfisr.lat(1,:,6),pfisr.lon(1,:,6));
plot_slice(pfisrNe.data,pfisr.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),pfisrNe.alt(1,:,:),'26 Mar 2008 11:00',122);
colorbar;


%%
figure;
scatter(pfisrNe.lat(1,:,86),pfisrNe.lon(1,:,86)); hold on; scatter(pfisr.lat(1,:,6),pfisr.lon(1,:,6));
plot_slice(log10(pfisr.energyFlux),pfisr.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),pfisr.zEnergyBin(1,:,:)/1000,'26 Mar 2008 11:00',3);
colorbar;
% caxis([64.99 65]);
hold on;


function [r,time,time2_corr] = calculate_correlation(data1,time1,data2,time2)
    
    time = time1; 
    nTime = length(time1);
    time2_corr = zeros(nTime,1);
    nZ = size(data1,3);
    if size(data1,3)~=size(data2,3)
        error('nZ values are not the same');
    end
    for iTime = 1:1:nTime
        iTime2 = find_time(time2,datestr(time(iTime)));
        for iZ = 1:1:nZ
            r(iTime,iZ) =  corr2(data1(iTime,:,iZ),data2(iTime2,:,iZ));
        end
        time2_corr(iTime) = time2(iTime2);
    end
    
end

function [r,time,time2_corr] = calculate_correlation_v2(data1,time1,data2,time2)
    
    time = time1; 
    nTime = length(time1);
    time2_corr = zeros(nTime,1);
    nZ = size(data1,3);
    if size(data1,3)~=size(data2,3)
        error('nZ values are not the same');
    end
    for iTime = 1:1:nTime
        data2_temp = interp1(time2,data2,time1(iTime),'linear');
        for iZ = 1:1:nZ
            r(iTime,iZ) =  corr2(data1(iTime,:,iZ),data2_temp(1,:,iZ));
        end
        time2_corr(iTime) = time(iTime);
    end
    
end


function generate_panels(h,data,time,lat,lon,zAxis,timeArrayStr,...
    zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr,xLimArr,yLimArr)

% zAxis should be 1/1000 of zValueBinArray
%  [eV]                       [keV]
%  [m]                        [km]

if nargin < 14
    yLimArr = [64.8 65.6];
end

if nargin < 13
    xLimArr = [-148.2 -146.8];
end

nTStr = length(timeArrayStr);
for i = 1:1:nTStr
    timeIndx(i) = find_time(time,timeArrayStr(i));
end
timeBinArray = time(timeIndx);

nZ = length(zValueBinArray);
nT = length(timeBinArray);

p = panel();

panelHeight = 20; %mm
panelBreadth = 20; %mm
rowHeights = repmat({{panelHeight}},1,nZ);
colBreadths = repmat({{panelBreadth}},1,nT);

p.pack(rowHeights,colBreadths);
p.de.margin = 4;
p.marginleft = 35;
p.margintop = 20;
title(p, titleStr);

for iTime = 1:1:nT
    for iZ = 1:1:nZ
        p(iZ,iTime).select();
        
        if iTime ~= 1
            set(gca,'YTickLabel',[]);
        else
            ylabel([num2str(zValueBinArray(iZ)),' ',zUnitStr]);
        end
        
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        else
            xlabel(datestr(timeBinArray(iTime),'HH:MM:ss'));
        end
        
           
        plot_slice(data,time,lat,lon,zAxis,...
            datestr(timeBinArray(iTime)),zValueBinArray(iZ));
        caxis(cLimArray(iZ,:));
        xlim(xLimArr);
        ylim(yLimArr);
        set(gca,'XTick',[-148,-147],'YTick',[65,65.5]);
        if iTime == nT
        colorbar_thin('Location','eastoutside','YLabel',colorUnitStr,'Width',1);
        end
    end
end
        
   
end

function generate_panels_asi(h,data,time,lat,lon,zAxis,timeArrayStr,...
    zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr,xLimArr,yLimArr)

if nargin < 14
    yLimArr = [-148.2 -146.8];
end

if nargin < 13
    xLimArr = [64.8 65.6];
end

% zAxis should be 1/1000 of zValueBinArray
%  [eV]                       [keV]
%  [m]                        [km]
nTStr = length(timeArrayStr);
for i = 1:1:nTStr
    timeIndx(i) = find_time(time,timeArrayStr(i));
end
timeBinArray = time(timeIndx);

nZ = length(zValueBinArray);
nT = length(timeBinArray);

p = panel();

panelHeight = 20; %mm
panelBreadth = 20; %mm
rowHeights = repmat({{panelHeight}},1,nZ);
colBreadths = repmat({{panelBreadth}},1,nT);

p.pack(rowHeights,colBreadths);
p.de.margin = 4;
p.marginleft = 35;

p.margintop = 20;
title(p, titleStr);

for iTime = 1:1:nT
    for iZ = 1:1:nZ
        p(iZ,iTime).select();
        
        if iTime ~= 1
            set(gca,'YTickLabel',[]);
        else
            ylabel([num2str(zValueBinArray(iZ)),' ',zUnitStr]);
        end
        
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        else
            xlabel(datestr(timeBinArray(iTime),'HH:MM:ss'));
        end
        
           
        plot_slice_asi(data,time,lat,lon,zAxis,...
            datestr(timeBinArray(iTime)),zValueBinArray(iZ));
        caxis(cLimArray(iZ,:));
        xlim(xLimArr);
        ylim(yLimArr);
        set(gca,'YTick',[-148,-147],'XTick',[65,65.5]);
        if iTime == nT
        colorbar_thin('Location','eastoutside','YLabel',colorUnitStr,'Width',1);
        end
    end
end
        
   
end

function plot_slice(data,time,lat,lon,zAxis,timeStr,zAxisValue)

timeIndx = find_time(time,timeStr);
dataSlice = data(timeIndx,:,:);
plot_2D_energy_slice_geodetic_v2019_v1(squeeze(dataSlice),...
    squeeze(lat),squeeze(lon),squeeze(zAxis),...
    [],zAxisValue,[],[],false,false,256);

end

function plot_slice_asi(data,time,lat,lon,zAxis,timeStr,zAxisValue)

timeIndx = find_time(time,timeStr);
dataSlice = data(timeIndx,:,:);
plot_2D_energy_slice_geodetic_v2019(dataSlice,...
    lat,lon,zAxis,...
    [],zAxisValue,[],[],false,false,256);

end

function plot_DASC_aer_v1 (data, az, el)
x = (90-el).*sind(360-az);
y = (90-el).*cosd(360-az);

h2=pcolor(x,y,data);
set(h2,'EdgeColor','none');
colorbar;
axis equal
% hold on; annotation('textarrow',[0.15 0.20],[0.8 0.8],'string','W');
% hold on; annotation('textarrow',[0.2 0.2],[0.85 0.90],'string','N');
hold on; text(-77,75,'N');
hold on; line([-75 -75],[70 60],'color','r');
hold on; text(-77,55,'S');

hold on; text(-69,65,'W');
hold on; line([-80 -70],[65 65],'color','r');
hold on; text(-85,65,'E');

end

function plot_scatter_aer(az,el,colorStr)
if nargin<3
    colorStr = 'r';
end
x = (90-el).*sind(360-az);
y = (90-el).*cosd(360-az);
scatter(x(:),y(:),3,colorStr);
end

function [i_neighbourhood,j_neighbourhood] = get_neighbouring_pixels(vi,vj,asiaz,asiel,beamWidth)
[I, J] = ndgrid(1:length(asiaz),1:length(asiel));
i_neighbourhood = zeros(size(vi,1),size(vi,2),64);
j_neighbourhood = zeros(size(vi,1),size(vi,2),64);
asizn = 90-asiel;
for i = 1:1:size(vi,1)
    for j= 1:1:size(vi,2)
        az0 = asiaz(vi(i,j),vj(i,j));
        zn0 = asizn(vi(i,j),vj(i,j));
        % translation of a circle in polar coordinates
        indx = asizn.^2 + zn0.^2 - 2.*asizn.*zn0.*cosd(asiaz-az0)<=beamWidth.*0.5;
        i_neighbourhood(i,j,1:length(I(indx))) = I(indx);
        j_neighbourhood(i,j,1:length(I(indx))) = J(indx);
    end
end

end

function asi_undersampled = get_undersampled_asi(asi,I,J)
% asi  is the all sky image [nTime x nPixels x nPixels]
% I are indices of interest [nBeams x nAlt or nEnergy x col number]
% J are indices of interest [nBeams x nAlt or nEnergy x row number]
% I(p,q,:) - the corresponding pixels in asi are to be averaged



% Get linear index numbers for a slice of ASI in time
nSize = size(squeeze(asi(1,:,:)));
nTime = size(asi,1);
nB = size(I,1);
nE = size(I,2);
asi_undersampled = zeros(nTime,nB,nE);

for iB = 1:1:nB
    for iE = 1:1:nE
    linearInd(iB,iE).arr = calculate_linear_index(nSize, I, J, iB, iE); 
    end
end

% Slicing ASI in time 
for iT = 1:1:nTime
    ASI = squeeze(asi(iT,:,:));
    
    for iB=1:1:nB
        for iE = 1:1:nE
            asi_undersampled(iT,iB,iE) = mean(ASI(linearInd(iB,iE).arr));
        end
    end
    
end

end

function linearInd = calculate_linear_index (matrixSize, i_n, j_n,iB,iE)

linearInd = squeeze(sub2ind(matrixSize,i_n(iB,iE,i_n(iB,iE,:)>0),j_n(iB,iE,j_n(iB,iE,:)>0)));

end

function [ax2, h2] = plot_2D_energy_slice_geodetic_v2019( diffEnergyFlux,...
    latitude, longitude, zEnergyBin, timeNumPFISR, ...
    thisEnergy, latWidth, lonWidth, setMapOn,...
    setTimeLabelOn, imageSize)
%plot_2D_energy_slice_geodetic.m Plot 2D differential energy flux slices 
%from 4-D PFISR data sets on lat, long map
%--------------------------------------------------------------------------
%Input
%-----
% data          : arranged beam-wise
% -> flux       : differential number flux [nE x nTime]
% -> energyFlux : differential energy flux [nE x nTime]
% -> chi2       : Reduced chi-2 of the maximum entropy regression [1 x nTime]
% -> qInvert    : Production rate from inverted energy flux [nh x nTime]
% -> maxIter    : Maximum number of iterations [1 x nTime]
% -> energyBin  : Energy bin values [nE x 1]
% -> time       : Time array [nTime x 1]
% -> alt        : Altitude array [nh x 1]
% -> A          : Production rate vs. number flux matrix [nh x nE]
% -> qInput     : Production rate derived from electron density [nh x nTime]
% amisrData     : 
% -> site       :latitude, longitude and altitude (PFISR location)
% -> magBeamNo  :the beam number/ID that points along the mag. field line 
% magcoords     : arranged non-beam-wise [nh x 3 x nBeams]
% energyBin     : Energy bin values [nE x 1]
% nBeams        : Total number of beams
% timeNo        : Time number of the energy slice to be plotted
% altitude      : Altitude of projection of the energy slice
% energy        : Energy in keV of the differential energy flux to be plotted
% setMapOn      : True => Map axis on
% setTimeLabelOn: True => The time and energy values are printed on the
%                 plot
%--------------------------------------------------------------------------
%Output
%------
% h2 - plot handle
%
%% 
%----------------------------------------------------------------------------
% Modified: 2nd Feb 2018 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
if nargin < 11
    imageSize = 64; %Hard coded? based on cal file
end
if nargin < 10
    setTimeLabelOn = true;
end
if nargin < 9
    setMapOn = true;
end
if nargin < 8
    lonWidth = 4;
end
if nargin < 7
    latWidth = 2;
end
%% Generating lat, lon, h, energy coordinates
% diffEnergyFlux(energy,beams) - at a time instant
% lat, lon, energyBin, diffenergyflux - for all data points
%% Generating data slice
F = scatteredInterpolant(latitude(:), longitude(:), zEnergyBin(:), diffEnergyFlux(:),'nearest','none');

latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);

Vq = F({latq,lonq,thisEnergy*1000});
Vq(Vq<=0)=nan;



if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
        'MapLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(latq,lonq,(Vq)); 
else
    h2=pcolor(latq,lonq,(Vq)); 
end


set(h2,'EdgeColor','none');

if setTimeLabelOn==true
    hold on;
    textm(latLim(2), lonLim(2)+lonWidth/10, ['PFISR: ', num2str(thisEnergy),' keV'],'color','r');
    textm(latLim(2)-latWidth/20, lonLim(2)+lonWidth/10,...
        [datestr(timeNumPFISR,'HH:MM:SS'),' UT'],'color', 'r');
    hold off;
end

end

function [ax2, h2] = plot_2D_energy_slice_geodetic_v2019_v1( diffEnergyFlux,...
    latitude, longitude, zEnergyBin, timeNumPFISR, ...
    thisEnergy, latWidth, lonWidth, setMapOn,...
    setTimeLabelOn, imageSize)
%plot_2D_energy_slice_geodetic.m Plot 2D differential energy flux slices 
%from 4-D PFISR data sets on lat, long map
%--------------------------------------------------------------------------
%Input
%-----
% data          : arranged beam-wise
% -> flux       : differential number flux [nE x nTime]
% -> energyFlux : differential energy flux [nE x nTime]
% -> chi2       : Reduced chi-2 of the maximum entropy regression [1 x nTime]
% -> qInvert    : Production rate from inverted energy flux [nh x nTime]
% -> maxIter    : Maximum number of iterations [1 x nTime]
% -> energyBin  : Energy bin values [nE x 1]
% -> time       : Time array [nTime x 1]
% -> alt        : Altitude array [nh x 1]
% -> A          : Production rate vs. number flux matrix [nh x nE]
% -> qInput     : Production rate derived from electron density [nh x nTime]
% amisrData     : 
% -> site       :latitude, longitude and altitude (PFISR location)
% -> magBeamNo  :the beam number/ID that points along the mag. field line 
% magcoords     : arranged non-beam-wise [nh x 3 x nBeams]
% energyBin     : Energy bin values [nE x 1]
% nBeams        : Total number of beams
% timeNo        : Time number of the energy slice to be plotted
% altitude      : Altitude of projection of the energy slice
% energy        : Energy in keV of the differential energy flux to be plotted
% setMapOn      : True => Map axis on
% setTimeLabelOn: True => The time and energy values are printed on the
%                 plot
%--------------------------------------------------------------------------
%Output
%------
% h2 - plot handle
%
%% 
%----------------------------------------------------------------------------
% Modified: 2nd Feb 2018 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
if nargin < 11
    imageSize = 64; %Hard coded? based on cal file
end
if nargin < 10
    setTimeLabelOn = true;
end
if nargin < 9
    setMapOn = true;
end
if nargin < 8
    lonWidth = 4;
end
if nargin < 7
    latWidth = 2;
end
%% Generating lat, lon, h, energy coordinates
% diffEnergyFlux(energy,beams) - at a time instant
% lat, lon, energyBin, diffenergyflux - for all data points
%% Generating data slice

iEnergy=find_altitude(zEnergyBin(1,:),thisEnergy);
latitude = latitude(:,iEnergy);
longitude = longitude(:,iEnergy);
diffEnergyFlux = diffEnergyFlux(:,iEnergy);

F = scatteredInterpolant( longitude(:), latitude(:), diffEnergyFlux(:),'nearest','none');

latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);

Vq = F({lonq,latq});
Vq(Vq<=0)=nan;



if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
        'MapLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
    h2=pcolorm(latq,lonq,(Vq)); 
else
    h2=pcolor(lonq,latq,(Vq)'); 
end


set(h2,'EdgeColor','none');

if setTimeLabelOn==true
    hold on;
    textm(latLim(2), lonLim(2)+lonWidth/10, ['PFISR: ', num2str(thisEnergy),' keV'],'color','r');
    textm(latLim(2)-latWidth/20, lonLim(2)+lonWidth/10,...
        [datestr(timeNumPFISR,'HH:MM:SS'),' UT'],'color', 'r');
    hold off;
end

end

