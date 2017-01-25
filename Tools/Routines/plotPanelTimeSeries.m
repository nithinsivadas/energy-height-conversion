%%  Routine to plot time-series panels

% % Plotting Wave Energy Figure
% dataChoice = [9,10,11,14,13];
% plot_preset_data_panels(dataChoice);
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.pdf' -pdf -nocrop
% save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.svg');
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.png' -r600 -png -nocrop
% 
% 
% % Plotting Substorm overview figure
dataChoice = [1,3,15,2,7,4,8];
plot_preset_data_panels(dataChoice);
if computer=='NITHIN-SURFACE'
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\SubstormOverview.pdf' -pdf -nocrop
    save('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\SubstormOverview.svg');
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\SubstormOverview.png' -r600 -png -nocrop
else
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft 12 Plots/SubstormOverview.pdf' -pdf -nocrop
    save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft 12 Plots/SubstormOverview.svg');
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft 12 Plots/SubstormOverview.png' -r600 -png -nocrop
end

% Plotting Inverted energy flux
dataChoice = [5,6,8,16,17,18,19];
plot_preset_data_panels(dataChoice);
if computer=='NITHIN-SURFACE'
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\InvertedEnergyFlux_with_Parallel_Potential.pdf' -pdf -nocrop
    save('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\InvertedEnergyFlux_with_Parallel_Potential.svg');
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\InvertedEnergyFlux_with_Parallel_Potential.png' -r600 -png -nocrop
else
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/InvertedEnergyFlux_with_Parallel_Potential.pdf' -pdf -nocrop
    % save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/InvertedEnergyFlux.svg');
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/InvertedEnergyFlux_with_Parallel_Potential.png' -r600 -png -nocrop
end
% % Plotting Parallel Potential Alone
% dataChoice = [16,18,19];
% plot_preset_data_panels(dataChoice);
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/Parallel_Potential.pdf' -pdf -nocrop
% % save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/InvertedEnergyFlux.svg');
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/Parallel_Potential.png' -r600 -png -nocrop

% Plotting Wave Energy
dataChoice = [9,10,11,14,13];
plot_preset_data_panels(dataChoice);
if computer=='NITHIN-SURFACE'
   export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\WaveEnergy.pdf' -pdf -nocrop
    save('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\WaveEnergy.svg');
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Paper 1\Plots\Draft 12 Plots\WaveEnergy.png' -r600 -png -nocrop
else
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.pdf' -pdf -nocrop
    save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.svg');
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/WaveEnergy.png' -r600 -png -nocrop
end

%% Test
% dataChoice = 3;
% plot_preset_data_panels(dataChoice);