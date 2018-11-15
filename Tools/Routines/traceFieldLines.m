%% Plot Field Lines on 26 Mar 2008
clear all;
%% Initializing
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 11:00','26 Mar 2008 12:00');

%%
tic
timeStr = '26-Mar-2008 11:30:00';
t = datetime(timeStr,'InputFormat','dd-MMM-yyyy HH:mm:ss');

GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));

thisTime = datenum(timeStr);
thisMaginput = interp1(timeMaginput',maginput,thisTime);
kp = round(thisMaginput(1)/10)+1;
PARMOD = zeros(10,1);
PARMOD(1) = kp;
PARMODT96 = zeros(10,1);
thisMaginput(2) = 4;
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
    
    [~,indx] = max(Kc{i});
    magEq{i} = [XX{i}(indx),YY{i}(indx),ZZ{i}(indx)];
    
end

%%
figure;
colormap(copper);
for i = 1:1:nLat
%     plot3(XX{i}',YY{i}',ZZ{i}');
    patch([XX{i} nan],[YY{i} nan],[ZZ{i} nan],[1./Kc{i} nan],'FaceColor','none','EdgeColor','interp');
    hold on;
    plot3(magEq{i}(1),magEq{i}(2),magEq{i}(3),'.r');
    text(magEq{i}(1),magEq{i}(2),magEq{i}(3)+((-1).^i)*0.2,...
        strcat(string(num2str(lat(i)')),'^0N'),'HorizontalAlignment','right');
end
cb = colorbar();

% set(gca,'colorscale','log')
cb.Ruler.Scale = 'log';
cb.Ruler.MinorTick = 'on';
% cb.Limits = [0.1 100];
caxis([0.01 1]);
xlabel('X_G_S_M [R_E]');
ylabel('Y_G_S_M [R_E]');
zlabel('Z_G_S_M [R_E]');
cb.Label.String = 'Radius of Curvature [R_E]';
title([datestr(thisTime),' Model:T96']);
view([0,0]);