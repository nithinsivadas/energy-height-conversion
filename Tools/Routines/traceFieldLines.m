%% Plot Field Lines on 26 Mar 2008
clear all;
%% Initializing
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
thmData = process_themis_data(datestr('26 Mar 2008 11:00'), [initialize_root_path,'LargeFiles',filesep],...
        'tha,thd,the','state');
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\lossConeFluxThmADE20080326.mat');
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 09:00','26 Mar 2008 12:00');

%%
timeArr = datenum('26-Mar-2008 11:00'):1/(60*24):datenum('26-Mar-2008 11:30');
tic
for i=1:1:length(timeArr)
timeStr = datestr(timeArr(i));
az = 180;
el = 0; %side
% el = 90; %top
fontsize = 14;
if el == 90
    viewStr = 'Top';
elseif el == 0
    viewStr = 'Side';
end

viewArr = [az,el];
storeImageDir = 'G:\My Drive\Research\Projects\Paper 2\Data\Figures\';
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

lat = [65,65.2,65.4,65.6];
nLat = length(lat);
alt = 110;
lon = -147.45;

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


% POES trajectory
poesTimeMinStr = '26-Mar-2008 11:27';
poesTimeMaxStr = '26-Mar-2008 11:33';
poesTimeIndx = find_time(poes.time,poesTimeMinStr):1:find_time(poes.time,poesTimeMaxStr);
poestime = poes.time(poesTimeIndx);
for i = 1:1:length(poestime)
    poesGSM(i,:) = onera_desp_lib_coord_trans([850,poes.lat(poesTimeIndx(i)),poes.lon(poesTimeIndx(i))],...
        [0 2],poestime(i));
end
poesColor = poes.mep0E2(poesTimeIndx);

% THM Coordinates

% THM trajectory
thmTimeMinStr = '26-Mar-2008 08:00';
thmTimeMaxStr = '26-Mar-2008 12:00';
energy = 100; %keV
thaTIndx = find_time(padData.tha.time,thmTimeMinStr):1:find_time(padData.tha.time,thmTimeMaxStr);
thaTime = padData.tha.time(thaTIndx);
thaGSM = interp1(thmData.tha.state.time,thmData.tha.state.XYZ_GSM,thaTime);
thaColor = interp1(padData.tha.energyBin,padData.tha.lcEfluxLi(thaTIndx,:)',energy*1000)';

thdTIndx = find_time(padData.thd.time,thmTimeMinStr):1:find_time(padData.thd.time,thmTimeMaxStr);
thdTime = padData.thd.time(thdTIndx,:);
thdGSM = interp1(thmData.thd.state.time,thmData.thd.state.XYZ_GSM,thdTime);
thdColor = interp1(padData.thd.energyBin,padData.thd.lcEfluxLi(thdTIndx,:)',energy*1000)';

theTIndx = find_time(padData.the.time,thmTimeMinStr):1:find_time(padData.the.time,thmTimeMaxStr);
theTime = padData.the.time(theTIndx);
theGSM = interp1(thmData.the.state.time,thmData.the.state.XYZ_GSM,theTime);
theColor = interp1(padData.the.energyBin,padData.the.lcEfluxLi(theTIndx,:)',energy*1000)';

% Magnetic Field Lines
% set(0,'Display
h=figure('visible','off');
ax1=axes;
ax1.Color = 'none';
colormap(ax1,copper);
for i = 1:1:nLat
%     plot3(XX{i}',YY{i}',ZZ{i}');
    patch([XX{i} nan],[YY{i} nan],[ZZ{i} nan],[1./Kc{i} nan],...
        'FaceColor','none','EdgeColor','interp','LineWidth',2);
    hold on;
    plot3(magEq{i}(1),magEq{i}(2),magEq{i}(3),'.r');
    text(magEq{i}(1),magEq{i}(2),magEq{i}(3)+((-1).^i)*0.2,...
        strcat(string(num2str(lat(i)')),'^0N'),'HorizontalAlignment','left');
    if strcmp(viewStr,'Side')
        text(magEq{i}(1),magEq{i}(2),2+((-1).^i)*0.2,...
            KEstr(i),'HorizontalAlignment','left','Color',[0.5 0.5 0.5]);
    end
end
cb = colorbar();
ax1.XLim = [-40 0];
ax1.ZLim = [-3 +3];
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
colormap(ax2,get_colormap('k','r'));
patch([poesGSM(:,1)' nan],[poesGSM(:,2)' nan],[poesGSM(:,3)' nan],[poesColor' nan],...
        'FaceColor','none','EdgeColor','interp','LineWidth',3);
cb2 = colorbar('Location','southoutside');
cb2.Position(2) = 0.2;
cb2.Position(3) = 0.25;
cb2.Position(4) = 0.01;
cb2.Label.String = 'NOAA17 EFlux >100 keV [Counts/s]';
ax2.CLim = [0 200];

% Plot THEMIS-A
    ax3 = axes;
    colormap(ax3,get_colormap('k','g'));
    patch([thaGSM(:,1)' nan],[thaGSM(:,2)' nan],[thaGSM(:,3)' nan],[thaColor' nan],...
            'FaceColor','none','EdgeColor','interp','LineWidth',3);
    hold on;
    thisTIndx = find_time(thaTime,'26-Mar-2008 08:00');
    time_label_3D(thaGSM(thisTIndx,1),thaGSM(thisTIndx,2),thaGSM(thisTIndx,3),thaTime(thisTIndx));
    thisTIndx = find_time(thaTime,'26-Mar-2008 12:00');
    time_label_3D(thaGSM(thisTIndx,1),thaGSM(thisTIndx,2),thaGSM(thisTIndx,3),thaTime(thisTIndx));
    thisTIndx = find_time(thaTime,datestr(thisTime));
    plot3(thaGSM(thisTIndx,1),thaGSM(thisTIndx,2),thaGSM(thisTIndx,3),'*m');

    cb3 = colorbar('Location','southoutside');
    cb3.Position = cb2.Position;
    cb3.Label.String = ['THEMIS-A Loss-cone EFlux ',num2str(energy),' keV'];
    cb3.Position(1) = 0.5;
    ax3.CLim = [0 1.4*10^4];

    % Plot THEMIS-D
    ax4 = axes;
    colormap(ax4,get_colormap('k','b'));
    patch([thdGSM(:,1)' nan],[thdGSM(:,2)' nan],[thdGSM(:,3)' nan],[thdColor' nan],...
        'FaceColor','none','EdgeColor','interp','LineWidth',3);
    hold on;
  
    thisTIndx = find_time(thdTime,thisTime);
    plot3(thdGSM(thisTIndx,1),thdGSM(thisTIndx,2),thdGSM(thisTIndx,3),'*m');
    
    
    cb4 = colorbar('Location','southoutside');
    cb4.Position = cb2.Position;
    cb4.Position(2) = 0.1;
    cb4.Label.String = ['THEMIS-D Loss-cone EFlux ',num2str(energy),' keV'];
    ax4.CLim = [0 5000];

    % Plot THEMIS-E
    ax5 = axes;
    colormap(ax5,get_colormap('k','c'));
    
    patch([theGSM(:,1)' nan],[theGSM(:,2)' nan],[theGSM(:,3)' nan],[theColor' nan],...
        'FaceColor','none','EdgeColor','interp','LineWidth',3);
    hold on;
  
    thisTIndx = find_time(theTime,thisTime);
    plot3(theGSM(thisTIndx,1),theGSM(thisTIndx,2),theGSM(thisTIndx,3),'*m');     

    cb5 = colorbar('Location','southoutside');
    cb5.Position = cb2.Position;
    cb5.Position(1) = 0.5;
    cb5.Position(2) = 0.1;
    cb5.Label.String = ['THEMIS-E Loss-cone EFlux ',num2str(energy),' keV'];
    ax5.CLim = [0 2500];

ax1.YLim =[min([ax3.YLim(1) ax4.YLim(1) ax5.YLim(1)]) max([ax3.YLim(2),ax4.YLim(2),ax5.YLim(2)])]; % Important
    
resize_figure(h,200,500);
ax1.OuterPosition = [0 0.25 1 0.75];
ax2 = copy_axes_properties(ax1,ax2);
ax3 = copy_axes_properties(ax1,ax3);
ax4 = copy_axes_properties(ax1,ax4);
ax5 = copy_axes_properties(ax1,ax5);

set(findall(gcf,'-property','FontSize'),'FontSize',fontsize);
export_fig(strcat(storeImageDir,imageName,'.png'),'-r300','-png','-nocrop');
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