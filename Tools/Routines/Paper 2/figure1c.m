%% For Paper 2, Figure 1c
% Generate figure with satellite tracks, and Themis All Sky Cameras
%% Load Data
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';


%% LOAD Digital MSP
dataMSP = get_msp_data('C:\Users\nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 2\Data\msp_vvel.mat','26 Mar 2008 10:00','26 Mar 2008 12:00');
pkrMLAT = 65.2639;
pkrMLON = -94.1995; 
h0=0.693;
elCutOff = 22.5;
el=fliplr(dataMSP.el+90);

%% Initalize

timeMinStr = '26 Mar 2008 11:20:18';
timeMaxStr = '26 Mar 2008 11:20:18';
time = datenum(timeMinStr):5/(24*60*60):datenum(timeMaxStr);
nTime = length(time);
latLim = [63,67];
lonLim = [-153,-141];
deltaLat = 1;
deltaLon = 5;
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\Figure1b\';
for i = 1:1:nTime
h=figure('visible','on');
% resize_figure(h,250,250);
[ax1,axm1]=combine_2D_plots_v2(fileStr,h,...
    'maps',{'OpticalImage','EnergyFluxMap'},...
    'sites',{'pokerFlat','pokerFlat'},...
    'thisTime',time(i),...
    'latLim',latLim,...
    'lonLim',lonLim,...
    'elCutOff',22.5,...
    'deltaLat',deltaLat,...
    'deltaLon',deltaLon,...
    'opticalLim',[250 400],...
    'peakIonizationAltitude',85,...
    'setStoreImage',false);

% Proton MSP Emissions
projMAlt = 110-h0; %km
[glat,glon,galt] = convert_mlatlon_to_glatlon(zeros(size(el)),el,projMAlt,pkrMLAT,pkrMLON,h0,time(i));

ax2=axes;
axm2 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax2.Color = 'none';
ax2.Visible = 'off';
cmap1 = colormap(ax2,get_colormap('w','m'));

thisTimeStr = datestr(time(i));
thisTimeIndxMSP = find_time(dataMSP.time,thisTimeStr);
scatterm(glat',glon'+0.7...
    ,7,dataMSP.intensity4861(:,thisTimeIndxMSP)','filled');

cb2 = colorbar('Location','southoutside');
cb2.Position(2) = 0.1;
cb2.Position(3) = 0.1;
cb2.Position(4) = 0.01;
cb2.Label.String = {'H_\beta 4861 A^0'};
ax2.CLim = [10 40];

% Blue emissions
projMAlt = 110-h0; %km
[glat,glon,galt] = convert_mlatlon_to_glatlon(zeros(size(el)),el,projMAlt,pkrMLAT,pkrMLON,h0,time(i));

ax3=axes;
axm3 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax3.Color = 'none';
ax3.Visible = 'off';
cmap2 = colormap(ax3,get_colormap('w','r'));

scatterm(glat',glon'+0.5,...
    7,dataMSP.intensity4278(:,thisTimeIndxMSP)'./(dataMSP.intensity5577(:,thisTimeIndxMSP)'),...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);

cb3 = colorbar('Location','southoutside');
cb3.Position(1) = 0.35;
cb3.Position(2) = 0.1;
cb3.Position(3) = 0.1;
cb3.Position(4) = 0.01;
cb3.Label.String = {'N_2^+/O^+'};
ax3.CLim = [0.1 0.3];

ax4=axes;
axm4 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax4.Color = 'none';
ax4.Visible = 'off';
cmap3 = colormap(ax4,get_colormap('w','b'));

scatterm(glat',glon'+0.5,...
    6,dataMSP.intensity4278(:,thisTimeIndxMSP)',...
    'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);

cb4 = colorbar('Location','southoutside');
cb4.Position(1) = 0.55;
cb4.Position(2) = 0.1;
cb4.Position(3) = 0.1;
cb4.Position(4) = 0.01;
cb4.Label.String = {'N_2^+ 4278 A^0'};
ax4.CLim = [0 200];

% Green emissions
projMAlt = 110-h0; %km
[glat,glon,galt] = convert_mlatlon_to_glatlon(zeros(size(el)),el,projMAlt,pkrMLAT,pkrMLON,h0,time(i));

ax5=axes;
axm5 = axesm('lambertstd','MapLatLimit',latLim,'MapLonLimit',lonLim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
    'PLineLocation',deltaLat,'MLineLocation',deltaLon);

ax5.Color = 'none';
ax5.Visible = 'off';
cmap2 = colormap(ax5,get_colormap('w','k'));

scatterm(glat',glon'+0.6,...
    7,dataMSP.intensity5577(:,thisTimeIndxMSP)',...
    'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);

cb5 = colorbar('Location','southoutside');
cb5.Position(1) = 0.75;
cb5.Position(2) = 0.1;
cb5.Position(3) = 0.1;
cb5.Position(4) = 0.01;
cb5.Label.String = {'O^+ 5577 A^0'};
ax5.CLim = [200 1000];

ax2 = copy_axes_properties(ax1,ax2);
ax3 = copy_axes_properties(ax1,ax3);
ax4 = copy_axes_properties(ax1,ax4);
ax5 = copy_axes_properties(ax1,ax5);


imageName = ['figure1b_',datestr(time(i),'HH_MM_SS')];
% export_fig(strcat(storeImageDir,imageName,'_with_POES.png'),'-r600','-png','-nocrop');
end






function plot_time_markers(time,lat,lon,markTimeArr,orientation,tickColor,textColor)
    if nargin<7
        textColor = 'k';
    end
    if nargin <6
        tickColor = 'k';
    end
        
    if strcmpi(orientation(1:5),'North')
%         dlat = +0.2;
        dlat = +0.05;
    elseif strcmpi(orientation(1:5),'South')
%         dlat = -0.2;
        dlat = -0.05;
    end
    if strcmpi(orientation(6:9),'East')
%         dlon = +0.5;
        dlon = +0.1;
    elseif strcmpi(orientation(6:9),'West')
%         dlon = -3;
        dlon = -1;
    end 
    
    for i = 1:length(markTimeArr)
        thisTimeIndx = find_time(time,datestr(markTimeArr(i)));
        plotm(lat(thisTimeIndx),lon(thisTimeIndx),'.','Color',tickColor);
        textm(lat(thisTimeIndx)+dlat,lon(thisTimeIndx)+dlon...
            ,datestr(markTimeArr(i),'HH:MM'),...
            'Color',textColor);
                %         create_perp_line(lat,lon,thisTimeIndx,orientation(1:5));
    end
    
end

function create_perp_line(lat,lon,markTimeArrIndx,orientation,tickLength)
    if nargin<6
        tickLength = 0.1;
    end
    
    if strcmpi(orientation(1:5),'North')
        dlat = 1;
    elseif strcmpi(orientation(1:5),'South')
        dlat = -1;
    end
    
    for i = length(markTimeArrIndx)
        thisTimeIndx = markTimeArrIndx(i);
        lat1 = lat(thisTimeIndx);
        lon1 = lon(thisTimeIndx);
        dxlon = mean(diff(lon(thisTimeIndx-1:thisTimeIndx+1)));
        dylat = gradient(lat(thisTimeIndx-1:thisTimeIndx+1),dxlon);
        slope = -1./dylat(2);
        c = lat1-slope*lon(1);
        lat2 = lat1+dlat*tickLength;
        lon2 = (lat2-c)./slope;
        hold on;
        plotm([lat1,lat2],[lon1,lon2],'k');
    end
    
end

function y = modSign(x)
    for i=1:1:length(x)
        if x(i)==0
            x(i) = +1;
        end
    end
    y = sign(x);
end

function [glat,glon,galt] = convert_mlatlon_to_glatlon(az,el,malt,mlat0,mlon0,h0,time0)
    C = define_universal_constants;
    slantRange = -C.RE.*sind(el) + modSign(el).*sqrt((C.RE^2).*(sind(el)).^2 + malt.*1000.*(malt.*1000+2.*C.RE)); 
%     slantRange = malt.*1000./sind(el);
%     slantRange(isinf(slantRange)) =slantRange(2);
% slant Range = RE*sin(el) +/- sqrt(RE^2 sin^2(el) + H(H+2RE); RE - radius
% of earth, H- altitude, el - elevation
    [mlat,mlon,h] = aer2geodetic(az,el,slantRange./1000,mlat0,mlon0,h0,wgs84Ellipsoid('km'));
%     mlat(el<22.5)=nan; mlat(el>180-22.5)=nan;
%     mlon(el<22.5)=nan; mlon(el>180-22.5)=nan;
%     h(el<22.5)=nan; h(el>180-22.5)=nan;
    
    tup=py.aacgmv2.wrapper.convert_latlon_arr(mlat,mlon,h-10.6,time0,'A2G');
    glat = double(py.array.array('d',py.numpy.nditer(tup{1})));
    glon = double(py.array.array('d',py.numpy.nditer(tup{2})));
    galt = double(py.array.array('d',py.numpy.nditer(tup{3})));
end