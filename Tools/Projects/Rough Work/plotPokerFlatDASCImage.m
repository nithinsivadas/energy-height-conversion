clear all;
close all;
azStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_AZ_10deg.FITS';
elStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_EL_10deg.FITS';
ASIDataStr = '/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_DASC_0000_20080326_114338.000.FITS';
% azStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_AZ_10deg.fits';
% elStr='/home/nithin/Documents/git-repos/Largefiles/Allsky_test/Cal_files_after_2011/PKR_DASC_20110112_EL_10deg.fits';
% ASIDataStr = '/home/nithin/Documents/git-repos/Largefiles/Allsky_test/PKR_DASC_0558_20151007_073824.971.FITS';

azOldRes=fitsread(azStr);
elOldRes=fitsread(elStr);

[dataNew, lat, lon, az_new, el_new, sensorloc, timeDASC] = DASC_aer_to_geodetic...
    (ASIDataStr, azOldRes, elOldRes,...
    512, 30, 110);

include_new_colormaps;

%%

figureHandle = figure;
resize_figure(figureHandle,210,148);
p=panel(figureHandle);
panelSize = 60; %in mm
p.pack({{90} {90}}, {{130}});
p.marginleft=5;
p.marginright=10;
p.fontsize=10;
p.de.margin=10;
p(1).margintop=10;
p(2).margintop = 10;
% p(2).margintop=6;

p(1,1).select();
plot_DASC_aer(dataNew, rotate_array(az_new,+18.13), el_new, 512); 
xlim([-110, 110]); ylim([-110,110]);
hold on; text(-100,100,'Geomagnetic North');
hold on; plot_grid_aer([0,90],[0, 30],'k');
hold on; plot_grid_aer(nan,[60, 90],'w');
colormap(viridis);
caxis([300 600]);

p(2,1).select();
plot_DASC_geodetic( dataNew, timeDASC, lat, lon, 512, [63 67], [-152 -143]);
colorbar;
caxis([300 600]);
colormap(viridis);
