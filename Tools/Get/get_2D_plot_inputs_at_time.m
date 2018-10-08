function plotData = get_2D_plot_inputs_at_time(inputH5FileStr, varargin)
% get_2D_plot_inputs_at_time.m Get inputs required to plot at a particular
% time instant for HDF5 file
%   Detailed explanation goes here
if nargin < 1
    error('inputH5FileStr not specified');
end

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
expectedMaps = {'OpticalImage','EnergyFluxMap','MagneticFieldMap','NoMap'};
expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};

addParameter(p,'plotModeStr','NoMap',@(x) any(validatestring(x,expectedMaps)));
addParameter(p,'magFieldModelStr','TS96',@(x) any(validatestring(x,expectedMagFieldModels)));
addParameter(p,'thisTimeIndx',1,validScalarPosNum);
addParameter(p,'plotData',struct());
addParameter(p,'getKc',false);
addParameter(p,'energySlice',100);

addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));

parse(p,inputH5FileStr,varargin{:});

plotData = p.Results.plotData;

switch p.Results.plotModeStr
    case 'OpticalImage'
        plotData.image = readh5_variable_at_time(inputH5FileStr,'ASI',...
            '/DASC/',p.Results.thisTimeIndx);
        plotData.thisTime = unix_to_matlab_time(readh5_variable_at_time(inputH5FileStr,'time',...
            '/DASC/',p.Results.thisTimeIndx));
    case 'EnergyFluxMap'
        plotData.diffEnergyFlux = readh5_variable_at_time(inputH5FileStr,...
            'energyFlux','/energyFluxFromMaxEnt/',p.Results.thisTimeIndx)';
        plotData.diffEnergyFlux(plotData.diffEnergyFlux<0) = 10^3;
        plotData.thisTime = readh5_variable_at_time(inputH5FileStr,...
            'time','/energyFluxFromMaxEnt/',p.Results.thisTimeIndx);
    case 'MagneticFieldMap'
        plotData.magEqCoordGEO = readh5_variable_at_time(inputH5FileStr,...
            'magEqCoordGEO',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx)';
        plotData.RE =(plotData.magEqCoordGEO(:,1).^2 + plotData.magEqCoordGEO(:,2).^2+...
                plotData.magEqCoordGEO(:,3).^2).^0.5;
        plotData.magFieldModelStr = p.Results.magFieldModelStr;
        try    
        plotData.Lm = abs(readh5_variable_at_time(inputH5FileStr,...
            'Lm',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx));
        plotData.Lstar = abs(readh5_variable_at_time(inputH5FileStr,...
            'Lstar',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx));
        catch ME;
        end
        plotData.thisTime = unix_to_matlab_time(readh5_variable_at_time(inputH5FileStr,...
            'time',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx));

        if p.Results.getKc
            plotData.BgeoEq = readh5_variable_at_time(inputH5FileStr,...
            'BgeoEq',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx)';
            plotData.BmagEq = readh5_variable_at_time(inputH5FileStr,...
            'BmagEq',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx);
            plotData.gradBmagEq = readh5_variable_at_time(inputH5FileStr,...
            'gradBmagEq',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx)';
            plotData.diffBEq = permute(readh5_variable_at_time(inputH5FileStr,...
            'diffBEq',['/magneticMap/',p.Results.magFieldModelStr,'/'],p.Results.thisTimeIndx),[3,2,1]);
            plotData.Kc = get_isotropic_boundary(plotData.BgeoEq,plotData.BmagEq,...
            plotData.gradBmagEq,plotData.diffBEq,p.Results.energySlice);
        end
        
    otherwise
        error('No maps or incorrect plotMode');
end



end

