%% Routine to test plotting combined 2D plots
%% Initializing
clear all;
inputH5FileStr='G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-energyFlux_5_Oct.h5';
magFieldModelStr = 'TS96';

%% Generating Inputs
pfisrData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','EnergyFluxMap','energySlice',100);

dascData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','OpticalImage');

magData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','MagneticFieldMap');
%% Combine more than one plot for a particular time instant
h=figure;
combine_2D_plots(inputH5FileStr,h,...
    'map1','OpticalImage','map3','EnergyFluxMap','map2','MagneticFieldMap',...
'map1Data',dascData,'map3Data',pfisrData,'map2Data',magData,...
'energySlice',100,'thisTime',datenum('26 Mar 2008 12:14'),'contourLineArray',1:30,...
'contourLabelArray',[5,7,9,12,15,20,25,30],...
'latLim',[63 67],'lonLim',[-153 -143],'opticalLim',[300 450],...
'plotContours','RE');
% 'figureLength',216,'figureBreadth',279);
