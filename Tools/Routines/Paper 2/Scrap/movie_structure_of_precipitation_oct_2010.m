%% Figure 1e (Variation of structure)
clear all;
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20101018.001_bc_15sec-full_v110.h5';
h5FileStrNe = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20101018.001_bc_15sec-full_v110.h5';
%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');

pfisrData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','EnergyFluxMap');
%%
dascImage = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
pfisrEnergyFluxImage = permute(h5read(h5FileStr,'/energyFluxFromMaxEnt/energyFlux'),[3 2 1]);
pfisrNeImage = permute(h5read(h5FileStr,'/inputData/Ne'),[3 2 1]);

%% Trim data
timeMinStr = '18 Oct 2010 07:30';
timeMaxStr = '18 Oct 2010 09:00';
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
% amisr = read_amisr(h5FileStrNe);
%%
% amisrData = aer_to_field_aligned_coords(amisr,60);
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

%% Global plotting variables

gl.timeArrayStr = string(['18 Oct 2010 07:45:00';...
                       '18 Oct 2010 08:00:00';...
                       '18 Oct 2010 08:08:18';...
                       '18 Oct 2010 08:12:31';...
                       '18 Oct 2010 08:18:05';...
                       '18 Oct 2010 08:20:55';...
                       '18 Oct 2010 08:30:00']);
gl.energyArray = [3,8,30,100];
gl.altitudeArray = [120,108,95,85];

%% Figure of Energy Flux 
h1 = figure;
resize_figure(h1,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.energyArray;
cLimArray = [7,9;...
            7,9;...
            7,9;...
            7,9];
zUnitStr = 'keV';
colorUnitStr = '[eV/m^2 sr s eV]';
titleStr = 'Energy flux from PFISR';

generate_panels(h1,log10(pfisr.energyFlux),pfisr.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
    pfisr.zEnergyBin(1,:,:)/1000,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
colormap('inferno');


%% Electron density (descrepancy between the coordinates of electron density and energy flux)
h2 = figure;
resize_figure(h2,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.altitudeArray;
cLimArray = [9.5,10.5;...
            9.5,10.5;...
            9.5,10.5;...
            9.5,10.5];
zUnitStr = 'km';
colorUnitStr = '[m^-^3]';
% uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1).*1000;
titleStr = 'Ne from PFISR';
pfisrNe.data(pfisrNe.data<=0) = nan;
generate_panels(h2,log10(pfisrNe.data),pfisr.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
    pfisrNe.alt(1,:,:),timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
colormap(magma);

%% Figure of asi undersampled/ energy flux
h3 = figure;
resize_figure(h3,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.energyArray;
cLimArray = [370,400;...
            370,400;...
            370,400;...
            370,400];
zUnitStr = 'keV';
colorUnitStr = '[a.u.]';
titleStr = 'ASI undersampled at PFISR energy flux coords';
generate_panels(h3,asi.eflux,asi.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
    pfisrData.zEnergyBin/1000,timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);
% colormap(get_colormap('k','g'));
colormap(viridis);

%% Figure of asi undersampled/ Ne
h4 = figure;
resize_figure(h4,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.altitudeArray;
cLimArray = [370,400;...
            370,400;...
            370,400;...
            370,400];
zUnitStr = 'km';
colorUnitStr = '[a.u.]';
titleStr = 'ASI undersampled at PFISR Ne coords';
% uniformAlt = repmat((pfisrNe.alt(1,13,:)),1,size(pfisrNe.data,2),1).*1000;
generate_panels(h4,asi.ne,asi.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
    pfisrNe.alt(1,:,:),timeArrayStr,zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr);

% colormap(get_colormap('k','g'));
colormap(viridis);

%% Figure of asi/energy flux with high resolution
h5=figure; 
resize_figure(h5,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.energyArray;
cLimArray = [370,400;...
            370,400;...
            370,400;...
            370,400];
zUnitStr = 'keV';
colorUnitStr = '[a.u.]';
titleStr = ' ASI in PFISR Field of View / Energy Flux';
generate_panels_asi(h5, asi.data, asi.time(:,1,1),asi.lat(1,:,:),asi.lon(1,:,:),...
    pfisr.lat(1,:,:),pfisr.lon(1,:,:),pfisr.zEnergyBin(1,1,:)./1000,timeArrayStr,...
    zValueBinArray, cLimArray, zUnitStr, colorUnitStr, titleStr);
% colormap(get_colormap('k','g'));
colormap(viridis);


%% Figure of asi/energy flux with high resolution
h6=figure; 
resize_figure(h6,125,225);
timeArrayStr = gl.timeArrayStr;
zValueBinArray = gl.altitudeArray;
cLimArray = [340,380;...
            340,380;...
            340,380;...
            340,380];
zUnitStr = 'km';
colorUnitStr = '[a.u.]';
titleStr = ' ASI in PFISR Field of View / Ne';
generate_panels_asi(h7, asi.data, asi.time(:,1,1),asi.lat(1,:,:),asi.lon(1,:,:),...
    pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),pfisrNe.alt(1,1,:),timeArrayStr,...
    zValueBinArray, cLimArray, zUnitStr, colorUnitStr, titleStr);
% colormap(get_colormap('k','g'));
colormap(viridis);

%% Find Correlation between ASI pixels and Ne/energy flux
[r_eflux,time_r_eflux, time2_r_eflux] = calculate_correlation_v2((pfisr.energyFlux),pfisr.time(:,1,1),(asi.eflux),asi.time(:,1,1));
%%
[r_ne,time_r_ne, time2_r_ne] = calculate_correlation_v2((pfisrNe.data),pfisr.time(:,1,1),(asi.ne),asi.time(:,1,1));
% [r_ne,time_r_ne, time2_r_ne] = calculate_correlation_v2((pfisrNe0.data),pfisrNe0.time(:,1,1),(asi.ne0),asi.time(:,1,1));
%% Find Correlation between Ne_@80km and the rest
altIndx = find_altitude(pfisrNe.alt(1,13,:),80);
pfisrne80km = repmat(interp_nans(pfisrNe.data(:,:,altIndx)')',1,1,size(pfisrNe.data,3));
[r_ne_ne,time_r_ne_ne] = calculate_correlation_v2(pfisrNe.data,pfisr.time(:,1,1),pfisrne80km,pfisr.time(:,1,1));

%%  Plot correlation between ASI pixels and Ne/ energy flux
totalPanelNo=2;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);
q=p(1);

q(1).select();
plot_2D_time_series(time_r_eflux,pfisrData.zEnergyBin(1,:)./1000,r_eflux',10/60,-1);
label_time_axis(time_r_eflux, false, 10/60);
set(gca,'YScale','log','YTick',[10,30,100,300],'YTickLabel',[10,30,100,300]);
ylabel('[keV]');
colormap(gca,get_colormap3('k','w','r'));
grid on;
colorbar_thin('YLabel',{'r_0(\phi_E & Optical'});
caxis([-1,1]);


q(2).select();
uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1);
plot_2D_time_series(time_r_ne,uniformAlt(1,13,:),r_ne',0.1,-1);
label_time_axis(time_r_ne, true, 10/60);
set(gca,'YScale','log','YTick',[70,80,100,120,200]);
colormap(gca,get_colormap3('k','w','r'));
% title('Correlation of Ne (mag-field-aligned) with ASI');
ylabel('[Km]');
colorbar_thin('YLabel',{'r_0 (N_e & Optical)'});
caxis([-1,1]);

%% Plot correlation between N_e@80 km and that at all other altitudes
totalPanelNo=1;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',25);
q=p(1);

q(1).select();
uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1);
plot_2D_time_series(time_r_ne_ne,uniformAlt(1,1,:),r_ne_ne',0.1,-1);
label_time_axis(time_r_ne_ne, true, 10/60);
set(gca,'YScale','log','YTick',[70,80,100,120,200]);
colormap(gca,get_colormap('k','w','r'));
% title('Correlation of Ne (mag-field-aligned) with ASI');
ylabel('[Km]');
colorbar_thin('YLabel',{'r_0 (N_e@80km & N_e)'});
caxis([-1,1]);
%% Testing if there is any time delay between the two images being correlated
figure;
plot(time_r_ne,(time2_r_ne(:)-time_r_ne(:))*24*3600);
label_time_axis(time_r_ne, true, 10/60);

%% Defining a polygon of radar field of view
figure;
beamlat = pfisr.lat(1,:,1);
beamlon = pfisr.lon(1,:,1);
scatter(beamlon',beamlat',15,1:26,'filled');
cmap = viridis(26);
b = num2str((1:26)');
hold 
text(beamlon,beamlat,cellstr(b),'Fontsize',10);
% cmap = flipud(cmap(1:10,:));
% cmap(1,:) = [1,1,1];
colormap(cmap);
colorbar;

beamIndx = [13,5,7,4,1,13];
xv = beamlon(beamIndx);
yv = beamlat(beamIndx);
x0 = mean(beamlon);
y0 = mean(beamlat);
[xv,yv] = push_point_away_from_center(xv,yv,x0,y0,0.1);
hold on;
plot(xv,yv);

%%
[in, on, pfisrFOV.lon, pfisrFOV.lat] = pfisr_field_of_view(asi.lat(1,:,:), asi.lon(1,:,:),...
            pfisrNe.lat(1,:,:), pfisrNe.lon(1,:,:),...
        pfisrNe.alt(1,1,:),85,0.1);
%%
correlation.r_eflux = r_eflux;
correlation.time_r_eflux = time_r_eflux;
correlation.r_ne = r_ne;
correlation.time_r_ne = time_r_ne;
%%
clc;
time = asi.time(:,1,1);
for iTime = 1:1:length(time)
% iTime = 100;
    100.*iTime./length(time)
    h=figure('visible','off');
    create_figure_sx(h,pfisr,pfisrNe,asi,correlation,pfisrFOV,gl,datestr(time(iTime)),...
        true,'G:\My Drive\Research\Projects\Paper 2\Videos\Oct18_2010\Movie\');
end
%%
create_video('G:\My Drive\Research\Projects\Paper 2\Videos\Oct18_2010\','Movie\','Movie_S1.avi');
%%
% clc;
% h=figure('visible','on');
% create_figure_sx(h,pfisr,pfisrNe,asi,correlation,pfisrFOV,gl,'26 Mar 2008 11:00',...
%         false,'G:\My Drive\Research\Projects\Paper 2\InternalReview1\Figure_edits\SliceMatrix\');

function h=create_figure_sx(h,pfisr,pfisrNe,asi,correlation,pfisrFOV,globalVar,timeStr,...
    setStoreImage,imageStoreDir)
resize_figure(h,215,325);
nZ = 4;
nT = 4;
p = panel();
p.pack({{80} {80}},{{150} {125}});

panelHeight = 20; %mm
panelBreadth = 20; %mm
rowHeights = repmat({{panelHeight}},1,nZ);
colBreadths = repmat({{panelBreadth}},1,nT);

p(1,2).pack(rowHeights,colBreadths);
p(1,2).de.margin = 4;
p(1,2).marginleft = 15;
p(1).de.margin = 4;
p.marginleft = 25;
p.margintop = 20;

p(2,1).pack({{30} {30}},{{125}});
p(2,1).de.margin = 4;
p(2,1,1).marginleft = 20;

p(2,2).pack({{50} },{{50} {50}});
p(2,2).marginleft = 20;
p(2,2).marginbottom = 10;
p(2,2).de.margin = 10;
p(2).margintop = 35;
p.select('all');


p(1,1).select();

time = asi.time(:,1,1);
iT = find_time(time,timeStr);
image = squeeze(asi.data(iT,:,:));
image = image.*2^-16;
avg = 367.*2^-16;
n=0.2;  
sigma = std2(image);
lowerLim = avg - n*sigma;
if lowerLim <=0
  lowerLim=0;
end
upperLim = avg + n*sigma;
if upperLim >=1
  upperLim=1;
end
image = imadjust(image, [lowerLim upperLim],[]);              
plot_DASC_geodetic(image, [], squeeze(asi.lat(iT,:,:)),squeeze(asi.lon(iT,:,:)),...
    1024,[63,67],[-153,-143],1,5);
colormap(get_colormap('k','g'));
caxis([0.3 1]);
colorbar_thin('YLabel','[a.u.]');

hold on;
plotm(pfisrFOV.lat,pfisrFOV.lon,'r');

% ASI - High resolution
zValueBinArray  = globalVar.altitudeArray;
zUnitStr = 'km';
cLimArray = [360,390;...
            360,390;...
            360,390;...
            360,390];
xLimArr = [-148.2 -146.8]; 
yLimArr = [64.8 65.5];
colorUnitStr = '[Counts]';
titleStr ={'ASI',colorUnitStr,datestr(asi.time(iT,1,1),'HH:MM:ss')};
for iZ = 1:1:nZ
    p(1,2,iZ).de.margin = 12;
    p(1,2,iZ,1).select();
    try
    [in, ~] = pfisr_field_of_view(asi.lat(1,:,:), asi.lon(1,:,:),...
        pfisrNe.lat(1,:,:), pfisrNe.lon(1,:,:),...
        squeeze(pfisrNe.alt(1,1,:)),zValueBinArray(iZ),0.05);
    plot_slice_asi(asi.data,asi.time(:,1,1),asi.lat(1,:,:),asi.lon(1,:,:),...
        timeStr,in);
    catch
    end
    caxis(cLimArray(iZ,:));
    xlim(xLimArr);
    ylim(yLimArr);
    set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'TickLength',[0.05 0.035]); 
    colorbar_thin('Location','eastoutside','Width',0.5);
    ylabel([num2str(zValueBinArray(iZ)),' ',zUnitStr]);
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        end
        if iZ == 1
            title(titleStr);
        end
end

% ASI - Low resolution
zValueBinArray  = globalVar.altitudeArray;
zUnitStr = 'km';
cLimArray = [360,390;...
            360,390;...
            360,390;...
            360,390];
colorUnitStr = '[Counts]';
titleStr ={'ASI-Undersampled',colorUnitStr,datestr(asi.time(iT,1,1),'HH:MM:ss')};
for iZ = 1:1:nZ
    p(1,2,iZ,2).select();
    try
    plot_slice(asi.ne,asi.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
        pfisrNe.alt(1,:,:),timeStr,zValueBinArray(iZ));
    catch
    end
    caxis(cLimArray(iZ,:));
    xlim(xLimArr);
    ylim(yLimArr);
    set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'YTickLabel',[],'TickLength',[0.05 0.035]); 
    colorbar_thin('Location','eastoutside','Width',0.5);
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        end
        if iZ == 1
            title(titleStr);
        end
end

% PFISR Ne
zValueBinArray  = globalVar.altitudeArray;
zUnitStr = 'km';
cLimArray = [9.5,10.5;...
            9.5,10.5;...
            9.5,10.5;...
            9.5,10.5];
colorUnitStr = '[m^-^3]';
iT1 = find_time(pfisr.time(:,1,1),timeStr);
titleStr ={'N_e',colorUnitStr,datestr(pfisr.time(iT1,1,1),'HH:MM:ss')};
for iZ = 1:1:nZ
    p(1,2,iZ,3).select();
%     iax = gca;
    plot_slice(real(log10(pfisrNe.data)),pfisr.time(:,1,1),pfisrNe.lat(1,:,:),pfisrNe.lon(1,:,:),...
        pfisrNe.alt(1,:,:),timeStr,zValueBinArray(iZ));
    caxis(cLimArray(iZ,:));
    xlim(xLimArr);
    ylim(yLimArr);
    set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'YTickLabel',[],'TickLength',[0.05 0.035]); 
    colormap(gca,viridis);
    colorbar_thin('Location','eastoutside','Width',0.5);
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        end
        if iZ == 1
            title(titleStr);
        end
end

% PFISR energyFLux
zValueBinArray  = globalVar.energyArray;
zUnitStr = 'keV';
cLimArray = [7,9;...
            7,9;...
            7,9;...
            7,9];
colorUnitStr = '[eV/m^2 sr s eV]';
titleStr ={'\phi(E)',colorUnitStr,datestr(pfisr.time(iT1,1,1),'HH:MM:ss')};
for iZ = 1:1:nZ
    p(1,2,iZ,4).marginleft = 15;
    p(1,2,iZ,4).select();
    plot_slice(real(log10(pfisr.energyFlux)),pfisr.time(:,1,1),pfisr.lat(1,:,:),pfisr.lon(1,:,:),...
        pfisr.zEnergyBin(1,:,:)/1000,timeStr,zValueBinArray(iZ));
    caxis(cLimArray(iZ,:));
    xlim(xLimArr);
    ylim(yLimArr);
    set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'YTickLabel',[],'TickLength',[0.05 0.035]); 
    colormap(gca,inferno);
    colorbar_thin('Location','eastoutside','YLabel',[num2str(zValueBinArray(iZ)),' ',zUnitStr],'Width',0.5);
        if iZ ~= nZ
            set(gca,'XTickLabel',[]);
        end
        if iZ == 1
            title(titleStr);
        end
end

% Correlation
timeIndx1 = find_time(correlation.time_r_eflux,timeStr);
timeIndx2 = find_time(correlation.time_r_ne,timeStr);

p(2,1,1,1).select();
plot_2D_time_series(correlation.time_r_eflux,pfisr.zEnergyBin(1,1,:)./1000,correlation.r_eflux',10/60,-1);
label_time_axis(correlation.time_r_eflux, false, 10/60);
set(gca,'YScale','log','YTick',[10,30,100,300],'YTickLabel',[10,30,100,300]);
ylabel('[keV]');
colormap(gca,get_colormap3('k','w','r'));
grid on;
colorbar_thin('YLabel',{'r_0 (\phi(E) & Optical)'});
caxis([-1,1]);
x = [correlation.time_r_eflux(timeIndx1) , correlation.time_r_eflux(timeIndx1)];
y = [min(pfisr.zEnergyBin(1,1,:)./1000), max(pfisr.zEnergyBin(1,1,:)./1000)];
hold on;
line(x,y,'Color','blue');

p(2,1,2,1).select();
uniformAlt = repmat(mean(pfisrNe.alt(1,:,:),2),1,size(pfisrNe.data,2),1);
plot_2D_time_series(correlation.time_r_ne,uniformAlt(1,1,:),correlation.r_ne',0.1,-1);
label_time_axis(correlation.time_r_ne, true, 10/60);
set(gca,'YScale','log','YTick',[70,80,100,120,200]);
colormap(gca,get_colormap3('k','w','r'));
% title('Correlation of Ne (mag-field-aligned) with ASI');
ylabel('[Km]');
colorbar_thin('YLabel',{'r_0 (N_e & Optical)'});
caxis([-1,1]);
x = [correlation.time_r_ne(timeIndx2) , correlation.time_r_ne(timeIndx2)];
y = [min(uniformAlt(1,1,:)), max(uniformAlt(1,1,:))];
hold on;
line(x,y,'Color','blue');

p(2,2,1,1).select();

plot(squeeze(pfisr.zEnergyBin(1,1,:)/1000),correlation.r_eflux(timeIndx1,:),'r');
set(gca,'XScale','log','XTick',[10,30,100,300],'XTickLabel',[10,30,100,300],...
    'TickLength',[0.05 0.035]);
xlim([1 1000]);
ylim([0 1]);
ylabel('r_0 (\phi(E) & Optical)');
xlabel('[keV]');
title(datestr(correlation.time_r_eflux(timeIndx1),'HH:MM:ss'));

p(2,2,1,2).select();

plot(squeeze(uniformAlt(1,1,:)),correlation.r_ne(timeIndx2,:),'k');
set(gca,'XScale','linear','XTick',[60,85,100,120,200],'XTickLabel',[60,85,100,120,200],...
    'TickLength',[0.05 0.035]);
ylim([0 1]);
xlim([60 200]);
ylabel('r_0 (N_e & Optical)');
xlabel('[km]');
title(datestr(correlation.time_r_ne(timeIndx2),'HH:MM:ss'));

    if setStoreImage == true
        export_fig(strcat(imageStoreDir,'Figure_',datestr(timeStr,'HH_MM_ss'),'.png'),...
            '-r300','-png','-nocrop');
        close(h);
    end

end

function [in,on,xv,yv] = pfisr_field_of_view(lat,lon,pfisrlat,pfisrlon,zAxis,zValue,dr)

zIndx = find_altitude(zAxis,zValue);
pfisrlat = pfisrlat(1,:,zIndx);
pfisrlon = pfisrlon(1,:,zIndx);

beamIndx = [13,5,7,4,1,13];
xv = pfisrlon(beamIndx);
yv = pfisrlat(beamIndx);
x0 = mean(pfisrlon);
y0 = mean(pfisrlat);
[xv,yv] = push_point_away_from_center(xv,yv,x0,y0,dr);
[in,on] = inpolygon(lon,lat,xv,yv);
end

function [x,y] = push_point_away_from_center(x,y,x0,y0,dr)
    x1 = x-x0;
    y1 = y-y0;
    
    [theta,rho] = cart2pol(x1,y1);
    [x,y] = pol2cart(theta,rho+dr);
    x = x+x0;
    y = y+y0;
    
end


% function [r,time,time2_corr] = calculate_correlation(data1,time1,data2,time2)
%     
%     time = time1; 
%     nTime = length(time1);
%     time2_corr = zeros(nTime,1);
%     nZ = size(data1,3);
%     if size(data1,3)~=size(data2,3)
%         error('nZ values are not the same');
%     end
%     for iTime = 1:1:nTime
%         iTime2 = find_time(time2,datestr(time(iTime)));
%         for iZ = 1:1:nZ
%             r(iTime,iZ) =  corr2(data1(iTime,:,iZ),data2(iTime2,:,iZ));
%         end
%         time2_corr(iTime) = time2(iTime2);
%     end
%     
% end

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
    yLimArr = [64.8 65.5];
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
        set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'TickLength',[0.05 0.035]);
        if iTime == nT
        colorbar_thin('Location','eastoutside','YLabel',colorUnitStr,'Width',1);
        end
    end
end
        
   
end

function generate_panels_asi(h,data,time,lat,lon,pfisrlat,pfisrlon,zAxis,timeArrayStr,...
    zValueBinArray,cLimArray,zUnitStr,colorUnitStr,titleStr,xLimArr,yLimArr)

if nargin < 16
    yLimArr = [64.8 65.5];
end

if nargin < 15
    xLimArr = [-148.2 -146.8];
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
        
        [in, ~] = pfisr_field_of_view(lat, lon,...
            pfisrlat(1,:,:), pfisrlon(1,:,:),...
        zAxis,zValueBinArray(iZ),0.05);
    
        plot_slice_asi(data,time,lat(1,:,:),lon(1,:,:),...
        datestr(timeBinArray(iTime)),in);   
        
        caxis(cLimArray(iZ,:));
        xlim(xLimArr);
        ylim(yLimArr);
        set(gca,'XTickMode','manual','YTickMode','manual',...
            'XTick',[-148,-147],'YTick',[65,65.5],'TickLength',[0.05 0.035]);
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

function plot_slice_asi(data,time,lat,lon,timeStr,in)

timeIndx = find_time(time,timeStr);
dataSlice = squeeze(data(timeIndx,:,:));
imagelat = squeeze(lat);
imagelon = squeeze(lon);
plot_DASC_geodetic_v1(dataSlice(squeeze(in))', [], imagelat(squeeze(in))', imagelon(squeeze(in))',...
    256,[],[],1,1);

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

function [ax1,h] = plot_DASC_geodetic_v1( dataNew, time, lat, lon,...
    imageSize, latLim, lonLim, deltaLat, deltaLon, label )
%% plot_DASC_geodetic Plot DASC optical image in geodetic coordinates on a map
%--------------------------------------------------------------------------
% Input
%------
% dataNew - 2-D optical data from DASC in geodetic coordinates [nCoordinates]
%           produced from DASC_aer_to_geodetic.m
% time    - Matlab time of the particular DASC image [nCoordinates]
% lat     - latitude coordinates [nCoordinates]
% lon     - longitude coordinates [nCoordinates]
% imageSize - the total pixel size of the output image e.g. 1024 
% latLim    - latitude limits e.g. [61 65]
% lonLim    - longitude limits e.g. [141.5 144.5]
%--------------------------------------------------------------------------
% Output
%-------
% ax1 - map axes
% h   - pcolor handle of color plot of the optical data
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
if nargin<8
    deltaLat = 1;
end

if nargin<9
    deltaLon = 5;
end

if nargin<10
    label = 'DASC';
elseif strcmp(label,'pokerFlat')
    label = 'PKFT DASC';
end

if imageSize == size(dataNew,1) 
    LAT = interp_nans(lat);
    LON = interp_nans(lon);
    Vq = dataNew; 
elseif imageSize == sqrt(length(dataNew))
    LAT = reshape(interp_nans(lat),imageSize,imageSize);
    LON = reshape(interp_nans(lon),imageSize,imageSize);
    Vq = reshape(dateNew,imageSize,imageSize);
else
    nanFlags = isnan(lat) | isnan(lon) | isnan(dataNew);
    lat(nanFlags) = [];
    lon(nanFlags) = [];
    dataNew(nanFlags) = [];
    % lat, lon, dataNew have to be column vectors!
    F = scatteredInterpolant(lat',lon',dataNew','natural','none'); %Modified on 3rd Oct 2018 - Transpose
    latq = linspace(min(lat),max(lat),imageSize);
    lonq = linspace(min(lon),max(lon),imageSize);
    [LAT, LON] = ndgrid(latq,lonq);
    Vq = F(LAT,LON);
end

% Plotting
% ax1=axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
%     'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
%     'PLineLocation',deltaLat,'MLineLocation',deltaLon);
% axis off
% load coastlines
% plotm(coastlat,coastlon,'Color',[0.5 0.5 0.5]);
% hold on;
h=pcolor(LON,LAT,(Vq)); 
set(h,'EdgeColor','none');
% if nargin>9
%     textm(latLim(2)+(latLim(2)-latLim(1))*0.05, lonLim(1) +(lonLim(2)-lonLim(1))*0.35, [char(upper(label)),': ',datestr(time,'HH:MM:SS'),' UT']);
% end
% textm(latLim(2)-0.19, lonLim(1)+0.1, ['DASC: ',datestr(time,'HH:MM:SS'),' UT']);
end



