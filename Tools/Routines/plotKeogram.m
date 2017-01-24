%% Plot Keogram of Optical, Electron density and Energy

clear all;
close all;
electronEnergy=100; %keV

azStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_AZ_10deg.FITS';
elStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_EL_10deg.FITS';
% ASIDataStr = '/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_DASC_0000_20080326_114338.000.FITS';
% azStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_AZ_10deg.fits';
% elStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_EL_10deg.fits';
% ASIDataStr = '/home/nithin/Documents/git-repos/Largefiles/Allsky_test/PKR_DASC_0558_20151007_073824.971.FITS';
MSPDataStr = '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Data/msp_vvel.mat';



azOldRes=fitsread(azStr);
elOldRes=fitsread(elStr);

% [dataNew, lat, lon, az_new, el_new, sensorloc, timeDASC] = DASC_aer_to_geodetic...
%     (ASIDataStr, azOldRes, elOldRes,...
%     512, 30, 110);
%% Prepare the Energy Spectra
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Data/2D_energy_spectra_26032008_800_to_1300Hr.mat');
[az,el,slant]=ned2aer(magcoords(:,2),magcoords(:,1),-magcoords(:,3));
zUp = magcoords(:,3);
aercoords = [az,el,slant,zUp];


include_new_colormaps;

%%

figureHandle = figure;

timeMinStr = '26 Mar 2008 11:00';
timeMaxStr = '26 Mar 2008 12:00';
dataMSP = extract_msp_data(MSPDataStr, timeMinStr, timeMaxStr);
%%
totalPanelNo=6;
p=panel(figureHandle);

p=panel();
p.pack(1);

panelSize = 30; %in mm
demargin = 4;
panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end;
p(1).pack(panelDefinition);

p.marginleft=35;
p.marginright=25;
p(1).de.margin=demargin;
% p.fontsize=12;
p.select('all');
timeTick=0.5;
  
resize_figure(figureHandle, panelSize*totalPanelNo+4*(totalPanelNo+1)+20);

dt=0.25;
 p(1,1).select();
 axPos1 = get(gca, 'position');
 plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 ylabel({['PFISR ',num2str(electronEnergy),' keV'],'Electrons','Elevation [deg]'});
 redmap=[1,1,1
     1,0.9,0.9
     1,0.8,0.8
     1,0.7,0.7
     1,0.6,0.6
     1,0.5,0.5
     1,0.4,0.4
     1,0.3,0.3
     1,0.2,0.2
     1,0.1,0.1
     1,0,0];
 colormap(gca,redmap);
 c=colorbar('eastoutside');
  ylabel(c,'log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]', 'FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos1);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
 set(c, 'Position',cPos);
  
 p(1,2).select();
 axPos = get(gca, 'position');
 plot_2D_time_series(dataMSP.time, dataMSP.el, log10(dataMSP.intensity4861), dt, 0, timeMinStr, timeMaxStr);
  set(gca,'YScale','linear');
 ylabel({'MSP 4861 A^0','H_\beta','Elevation [deg]'});
 colormap(gca,'viridis');
 c=colorbar('eastoutside');
 ylabel(c,'log_1_0 Intensity [a.u.]','FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
%  label_time_axis(dataMSP.time, true, dt, timeMinStr, timeMaxStr);
 % Overlaying 100 keV Electrons
 set(c, 'Position',cPos);
 axOverlay=axes('Position',get(gca,'Position'));
 h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 alpha color;
 colormap(gca,redmap);
 set(gca,'visible','off')
 
  p(1,3).select();
 axPos = get(gca, 'position');
 plot_2D_time_series(dataMSP.time, dataMSP.el, log10(dataMSP.intensity5577), dt, 0, timeMinStr, timeMaxStr);
  set(gca,'YScale','linear');
 ylabel({'MSP 5577 A^0','O I','Elevation [deg]'});
 colormap(gca,'viridis');
 c=colorbar('eastoutside');
 ylabel(c,'log_1_0 Intensity [a.u.]','FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
%  label_time_axis(dataMSP.time, true, dt, timeMinStr, timeMaxStr);
 % Overlaying 100 keV Electrons
 set(c, 'Position',cPos);
 axOverlay=axes('Position',get(gca,'Position'));
 h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 alpha color;
 colormap(gca,redmap);
 set(gca,'visible','off')
 
  p(1,4).select();
 axPos = get(gca, 'position');
 plot_2D_time_series(dataMSP.time, dataMSP.el, log10(dataMSP.intensity4278), dt, 0, timeMinStr, timeMaxStr);
  set(gca,'YScale','linear');
 ylabel({'MSP 4278 A^0','N_2^+','Elevation [deg]'});
 colormap(gca,'viridis');
 c=colorbar('eastoutside');
  ylabel(c,'log_1_0 Intensity [a.u.]','FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
%  label_time_axis(dataMSP.time, true, dt, timeMinStr, timeMaxStr);
 % Overlaying 100 keV Electrons
 set(c, 'Position',cPos);
 axOverlay=axes('Position',get(gca,'Position'));
 h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 alpha color;
 colormap(gca,redmap);
 set(gca,'visible','off')
 
  p(1,5).select();
 axPos = get(gca, 'position');
 plot_2D_time_series(dataMSP.time, dataMSP.el, log10(dataMSP.intensity6300), dt, 0, timeMinStr, timeMaxStr);
  set(gca,'YScale','linear');
 ylabel({'MSP 6300 A^0','O I','Elevation [deg]'});
 colormap(gca,'viridis');
 c=colorbar('eastoutside');
  ylabel(c,'log_1_0 Intensity [a.u.]','FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
 label_time_axis(dataMSP.time, true, dt, timeMinStr, timeMaxStr);
 % Overlaying 100 keV Electrons
 set(c, 'Position',cPos);
 axOverlay=axes('Position',get(gca,'Position'));
 h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 alpha color;
 colormap(gca,redmap);
 set(gca,'visible','off')
 
   p(1,6).select();
 axPos = get(gca, 'position');
 plot_2D_time_series(dataMSP.time, dataMSP.el, (dataMSP.intensity4278./dataMSP.intensity5577), dt, 0, timeMinStr, timeMaxStr);
  set(gca,'YScale','linear');
 ylabel({'Ratio 4278/5577 ','Hardness of Spectra','Elevation [deg]'});
 colormap(gca,'viridis');
 caxis([0 2]); 
 c=colorbar('eastoutside');
  ylabel(c,'Ratio','FontSize',8);
 % reducing color bar thickness
 set(gca,'position',axPos);
 cPos=get(c,'Position');
 cPos(3)=0.2*cPos(3);
 label_time_axis(dataMSP.time, true, dt, timeMinStr, timeMaxStr);
 % Overlaying 100 keV Electrons
 set(c, 'Position',cPos);
 axOverlay=axes('Position',get(gca,'Position'));
 h=plot_energyflux_keogram( data, aercoords, energyBin, nBeams, timeMinStr, timeMaxStr, 110, electronEnergy );
 alpha color;
 colormap(gca,redmap);
 set(gca,'visible','off')
 


