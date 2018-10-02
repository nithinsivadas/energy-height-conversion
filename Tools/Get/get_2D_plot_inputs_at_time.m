function plotData = get_2D_plot_inputs_at_time(inputH5FileStr, plotModeStr, thisTimeIndx, magFieldModelStr)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    error('inputH5FileStr not specified');
end

if nargin < 2 || isempty(plotModeStr)
    plotModeStr = 'Default';
end
   

switch plotModeStr
    case 'Optical Image'
        plotData.image = readh5_variable_at_time(inputH5FileStr,'ASI',...
            '/DASC/',thisTimeIndx);
        plotData.thisTime = unix_to_matlab_time(readh5_variable_at_time(inputH5FileStr,'time',...
            '/DASC/',thisTimeIndx));
    case 'Energy Flux Map'
        plotData.diffEnergyFlux = readh5_variable_at_time(inputH5FileStr,...
            'energyFlux','/energyFluxFromMaxEnt/',thisTimeIndx)';
        plotData.diffEnergyFlux(plotData.diffEnergyFlux<0) = 10^3;
        plotData.thisTime = unix_to_matlab_time(readh5_variable_at_time(inputH5FileStr,...
            'time','/DASC/',thisTimeIndx));
    case 'Magnetic Field Map'
        if nargin<4 || isempty(magFieldModelStr)
            magFieldModelStr = 'TS96';
        end
        plotData.magEqCoordGEO = readh5_variable_at_time(inputH5FileStr,...
            'magEqCoordGEO',['/magneticMap/',magFieldModelStr','/'],thisTimeIndx);
        plotData.thisTime = readh5_variable_at_time(inputH5FileStr,...
            'time',['/magneticMap/',magFieldModelStr','/'],thisTimeIndx);
    otherwise
        error('No or incorrect plotMode');
end



end

