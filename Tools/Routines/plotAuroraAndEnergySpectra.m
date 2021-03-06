clear all;
close all;
rootpath=initialize_root_path;
load ([rootpath,'energy-height-conversion/Tools/Projects/Paper 2/Data/2D_energy_spectra_26032008_800_to_1300Hr.mat']);
azStr=[rootpath,'Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_AZ_10deg.FITS'];
elStr=[rootpath,'Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_EL_10deg.FITS'];
include_new_colormaps;

ASIDataStr = [rootpath,'Largefiles/PokerFlat_DASC_08_03_26/PKR_DASC_0000_20080326_112558.000.FITS'];
% azStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_AZ_10deg.fits';
% elStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_EL_10deg.fits';
% ASIDataStr = '/home/nithin/Documents/git-repos/Largefiles/Allsky_test/PKR_DASC_0558_20151007_073824.971.FITS';

azOldRes=fitsread(azStr);
elOldRes=fitsread(elStr);

[dataNew, lat, lon, az_new, el_new, sensorloc, timeDASC] = DASC_aer_to_geodetic...
    (ASIDataStr, azOldRes, elOldRes,...
    512, 30, 110);



%% Create two axes
figureHandle = figure;
resize_figure(figureHandle,148,210); %A5 Paper Size
[axesHandleOptical, h1]=plot_DASC_geodetic(dataNew, timeDASC, lat, lon, 512, [63 67], [-152 -143]);

axesHandleEnergy = axes;
axesm('lambertstd','MapLatLimit',getm(axesHandleOptical,'MapLatLimit'),...
        'MapLonLimit',getm(axesHandleOptical,'MapLonLimit'),...
    'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
    'PLineVisible','off','MLineVisible','off')
axis off
timePFISRNo = find_time(data(1).time,datestr(timeDASC));
[h2]=plot_2D_energy_slice_geodetic(data, magcoords, energyBin, nBeams, timePFISRNo, 110, 30, false);

colormap(axesHandleOptical,'viridis');
colormap(axesHandleEnergy,'inferno');
cb1 = colorbar(axesHandleOptical,'eastoutside');
cb2 = colorbar(axesHandleEnergy,'westoutside');
ylabel(cb1, '[Rayleigh]');                  
ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
caxis(axesHandleOptical,[300 600]);
caxis(axesHandleEnergy,[8 10]);
% caxis(axesHandleEnergy,'auto');
alpha(axesHandleEnergy,0.5);
set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
linkaxes([axesHandleOptical,axesHandleEnergy]);

% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample Images/Sample_aurora_energy.pdf' -pdf -nocrop
% save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample_aurora_energy.svg');
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample Images/Sample_aurora_energy.png' -r600 -png -nocrop

