%% Event 2nd Feb 2019
clear all;
%%
fileStr = 'G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\PFISR\20190202.008_ac_1min-fitcal.h5';
% fileStr = 'G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\PFISR\20190302.001_bc_1min-fitcal.h5';
%%
load('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\19061_19033\19061_19033\19033\alt_Z.mat');
load('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\19061_19033\19061_19033\19033\Na_Z.mat');
load('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\19061_19033\19061_19033\19033\time.mat');

lidar.alt = alt_z;
lidar.nNa = n_na_u;
lidar.time = times_av;

%% beam1 txt
% lidar.alt=dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den_alt1.txt');
% times_av2=dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den_time1.txt');
% lidar.nNa = dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den1.txt');
%% beam2 txt
% lidar.alt=dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den_alt2.txt');
% times_av2=dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den_time2.txt');
% lidar.nNa = dlmread('G:\My Drive\Research\Projects\Collaborations\Katrina Bossert\LIDAR\02_Feb_2019\19033_den2.txt');

%%
data = read_amisr(fileStr);
%% Plot
% iB =2;
iB = 7;
h= figure;
timeMinStr = '02 Feb 2019 08:00';
timeMaxStr = '02 Feb 2019 13:00';
lidar.time = floor(datenum(timeMinStr))+times_av/24;
% resize_figure(h,
p = create_panels(h,'totalPanelNo',2,'demargin',4,'panelHeight',60,'panelBreadth',300);
p(1,1).select();
electronDensity = data.electronDensity;
electronDensity(electronDensity<=0) = 1;
plot_2D_time_series(data.time(1,:),data.altitude(:,iB),...
    log10(squeeze(electronDensity(:,iB,:))),1,0,timeMinStr, timeMaxStr);
colormap(inferno);
c=colorbar_thin();
ylabel(c,['log_1_0(N_e) Zenith Beam (No.',num2str(iB),')']);
ylabel({'PFISR N_e Densities','Altitude [km]'});
% label_time_axis(data.time(1,:),false,1,timeMinStr,timeMaxStr);
ylim([80 105]);
caxis([10 12]);
grid on;
title(datestr(floor(datenum(timeMinStr))));

p(1,2).select();
nNa = lidar.nNa;
nNa(nNa<=0) = 1;
plot_2D_time_series(lidar.time(1,:),lidar.alt(1,:),...
    (nNa),1,0,timeMinStr, timeMaxStr);
colormap(gca,viridis);
c=colorbar_thin();
ylabel(c,'Zenith Beam');
% ylabel(c,'Zenith Beam');
ylabel({'LIDAR Na Densities','Altitude [km]'});
label_time_axis(lidar.time(1,:),true,1,timeMinStr,timeMaxStr);
ylim([80 105]);
caxis([1e9 10e9]);