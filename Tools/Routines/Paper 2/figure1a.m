%% For Paper 2, Figure 1a
%% Plot Field Lines on 26 Mar 2008
clear all;
%% Initializing
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
thmData = process_themis_data(datestr('26 Mar 2008 11:00'), [initialize_root_path,'LargeFiles',filesep],...
        'tha,thd,the','state');
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');

%%
% POES trajectory
poesTimeMinStr = '26-Mar-2008 11:27';
poesTimeMaxStr = '26-Mar-2008 11:33';
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
end

% THM Coordinates

% THM trajectory
thmTimeMinStr = '26-Mar-2008 08:00';
thmTimeMaxStr = '26-Mar-2008 12:00';

thdTIndx = find_time(padData.thd.time,thmTimeMinStr):1:find_time(padData.thd.time,thmTimeMaxStr);
thdTime = padData.thd.time(thdTIndx,:);
thdGSM = interp1(thmData.thd.state.time,thmData.thd.state.XYZ_GSM,thdTime);

theTIndx = find_time(padData.the.time,thmTimeMinStr):1:find_time(padData.the.time,thmTimeMaxStr);
theTime = padData.the.time(theTIndx);
theGSM = interp1(thmData.the.state.time,thmData.the.state.XYZ_GSM,theTime);

%%
% timeArr = datenum('26-Mar-2008 11:00'):1/(60*24):datenum('26-Mar-2008 11:30');
timeArr = datenum('26-Mar-2008 11:19');
lat = [60,65.14,65.18,65.20,65.21,65.3];
setLabel =false;
fontsize = 14;
ZLim = [-3 +3];
XLim = [-15 0];

tic
for i=1:1:length(timeArr)
timeStr = datestr(timeArr(i));
az = 180;
el = 0; %side
% el = 90; %top

if el == 90
    viewStr = 'Top';
elseif el == 0
    viewStr = 'Side';
end

viewArr = [az,el];
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\Draft\';
imageName = ['TS96_',viewStr,datestr(datenum(timeStr),'HH_MM')];

t = datetime(timeStr,'InputFormat','dd-MMM-yyyy HH:mm:ss');

GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));
global GEOPACK1;
C = define_universal_constants;

thisTime = datenum(timeStr);
thisMaginput = interp1(timeMaginput',maginput,thisTime);
kp = round(thisMaginput(1)/10)+1;
PARMOD = zeros(10,1);
PARMOD(1) = kp;
PARMODT96 = zeros(10,1);
PARMODT96(1) = thisMaginput(5);
PARMODT96(2) = thisMaginput(2);
PARMODT96(3) = thisMaginput(6);
PARMODT96(4) = thisMaginput(7);

nLat = length(lat);
alt = 110;
lon = -147.47;

for i = 1:1:nLat

    thisLat = lat(i);
    xGEO = onera_desp_lib_rotate([alt, thisLat, lon],'gdz2geo');
    xGSM = onera_desp_lib_rotate(xGEO,'geo2gsm',thisTime);

    [~,~,~,XX{i},YY{i},ZZ{i},~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
        1.0,50,(6371.2+110)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');
    
    Kc{i} = geopack_find_curvature(XX{i},YY{i},ZZ{i});
    
    
    [maxKc,indx] = max(Kc{i});
    minRc = 1./maxKc;
    magEq{i} = [XX{i}(indx),YY{i}(indx),ZZ{i}(indx)];
    [BX,BY,BZ] = T96(0,PARMODT96,GEOPACK1.PSI,magEq{i}(1),magEq{i}(2),magEq{i}(3));
    KE(i) = ((((minRc.*C.RE ).*(C.e).*(BZ*10^-9)).^2).*(2^-7).*(C.me).^-1).*(10^-3).*(C.e).^-1; %keV
    if KE(i) < 1
        KEstr(i) = string([num2str(KE(i),'%3.2f'),' keV']);
    elseif KE(i) >= 1000
        KEstr(i) = string([num2str(KE(i)/1000,'%3.0f'),' MeV']);
    else
        KEstr(i) = string([num2str(KE(i),'%3.0f'),' keV']);
    end
end

% Magnetic Field Lines
h=figure('visible','on');
ax1=axes;
ax1.Color = 'none';
colormap(ax1,copper);
for i = 1:1:nLat
    patch([XX{i} nan],[YY{i} nan],[ZZ{i} nan],[1./Kc{i} nan],...
        'FaceColor','none','EdgeColor','interp','LineWidth',2);
    hold on;
    plot3(magEq{i}(1),magEq{i}(2),magEq{i}(3),'.r');
    if setLabel
        text(magEq{i}(1),magEq{i}(2),magEq{i}(3)+((-1).^i)*0.2,...
            strcat(string(num2str(lat(i)')),'^0N'),'HorizontalAlignment','left');
        if strcmp(viewStr,'Side')
            text(magEq{i}(1),magEq{i}(2),2+((-1).^i)*0.2,...
                KEstr(i),'HorizontalAlignment','left','Color',[0.5 0.5 0.5]);
        end
    end
end

    cb = colorbar();
    ax1.XLim = XLim;
    ax1.ZLim = ZLim;
    cb.Ruler.Scale = 'log';
    cb.Ruler.MinorTick = 'on';
    caxis([0.01 1]);
    xlabel('X_G_S_M [R_E]');
    ylabel('Y_G_S_M [R_E]');
    zlabel('Z_G_S_M [R_E]');
    cb.Label.String = 'Radius of Curvature [R_E]';
    title([datestr(thisTime),' Model:T96']);
    view(viewArr);


    % Plot POES
    
    hold on;
    ax2 = axes;
    plot3(poesGSM(:,1),poesGSM(:,2),poesGSM(:,3),'Color','r','LineWidth',2);
    ax2.Color = 'none';
    hold on;

    % Plot THEMIS-D
    ax4 = axes;
    plot3(thdGSM(:,1),thdGSM(:,2),thdGSM(:,3),'Color','b','LineWidth',2);
    
    hold on;
    thisTIndx = find_time(thdTime,thisTime);
    plot3(thdGSM(thisTIndx,1),thdGSM(thisTIndx,2),thdGSM(thisTIndx,3),'*m');
    ax4.Color = 'none';
    
    
    % Plot THEMIS-E
    ax5 = axes;
    plot3(theGSM(:,1),theGSM(:,2),theGSM(:,3),'Color','c','LineWidth',2);
    hold on;
    ax5.Color = 'none';
    thisTIndx = find_time(theTime,thisTime);
    plot3(theGSM(thisTIndx,1),theGSM(thisTIndx,2),theGSM(thisTIndx,3),'*m');     
 
    % Setting axes to scale
    ax1.YLim =[min([ax1.YLim(1),ax4.YLim(1),ax5.YLim(1)]) max([ax1.YLim(2),ax4.YLim(2),ax5.YLim(2)])]; % Important
    
    resize_figure(h,200,500);
    set(findall(gcf,'-property','FontSize'),'FontSize',fontsize);
    ax1.OuterPosition = [0 0.25 1 0.75];
    ax2 = copy_axes_properties(ax1,ax2);
    ax4 = copy_axes_properties(ax1,ax4);
    ax5 = copy_axes_properties(ax1,ax5);
    set(ax2,'visible','off')
    set(ax4,'visible','off')
    set(ax5,'visible','off')
    set(h, 'PaperPositionMode', 'auto');

export_fig(strcat(storeImageDir,imageName,'final_no_label.png'),'-r300','-png','-nocrop');
end
toc

function time_label_3D(x,y,z,time,format,color,verticalAlignment)

if nargin<6
    color = [0.5,0.5,0.5];
end
if nargin<5
    format = 'HH';
end

    text(x,y,z,{'_|',datestr(time,format)},'HorizontalAlignment','center','Color',color);
end