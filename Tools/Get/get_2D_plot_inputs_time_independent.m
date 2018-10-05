function plotData = get_2D_plot_inputs_time_independent(inputH5FileStr,varargin)
% get_2D_plot_inputs_at_time.m Get inputs required to plot at a particular
% time instant for HDF5 file
%   Detailed explanation goes here

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

expectedMaps = {'OpticalImage','EnergyFluxMap','MagneticFieldMap','NoMap'};
expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};

addParameter(p,'plotModeStr','NoMap',@(x) any(validatestring(x,expectedMaps)));
addParameter(p,'magFieldModelStr','TS96',@(x) any(validatestring(x,expectedMagFieldModels)));
addParameter(p,'timeNeutralAtmosphere',9999,validScalarPosNum);
addParameter(p,'energySlice',100,validScalarPosNum); %In keV

addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));

parse(p,inputH5FileStr,varargin{:});


if nargin < 1
    error('inputH5FileStr not specified');
end


switch p.Results.plotModeStr
    case 'OpticalImage'
        plotData.latitude = readh5_variable_at_time(inputH5FileStr,'lat',...
            '/DASC/',[]);
        plotData.longitude = readh5_variable_at_time(inputH5FileStr,'lon',...
            '/DASC/',[]);
        plotData.time = unix_to_matlab_time(h5read(inputH5FileStr,'/DASC/time'))';
    case 'EnergyFluxMap'
        
        magcoords = permute(readh5_variable_at_time(inputH5FileStr,...
            'magGeodeticLatLonAlt','/magneticFieldAlignedCoordinates/',[]),[3 2 1]);
              
        % Interpolate beam-wise at projectionAlt, and get the lat, lon values for [nBeamsxnEnergy]
        tempLatitude = squeeze(magcoords(1,:,:)); 
        tempLongitude = squeeze(magcoords(2,:,:));
        tempAltitude = squeeze(magcoords(3,:,:)); 
        nBeams = size(tempLatitude,1);
        plotData.zEnergyBin = repmat(readh5_variable_at_time(inputH5FileStr,...
            'energyBin','/energyFluxFromMaxEnt/',[])',nBeams,1);
        
        plotData.time = h5read(inputH5FileStr,'/energyFluxFromMaxEnt/time')';
        if p.Results.timeNeutralAtmosphere == 9999
            timeNeutralAtmosphere=median(plotData.time);
        else
            timeNeutralAtmosphere=p.Resilts.timeNeutralAtmosphere;
        end
        peakIonizationAlt = calculate_peak_altitude_of_ionization(p.Results.energySlice*1000,...
            timeNeutralAtmosphere,tempLatitude(1,:)',tempLongitude(1,:)',...
            tempAltitude(1,:)');
        fprintf(['PFISR Energy Flux Map projected at ',num2str(peakIonizationAlt),' km\n']);
        nBeams = size(tempLatitude,1);
        for iBeam = 1:1:nBeams
            pfisrLatitude(iBeam,1) = interp1(tempAltitude(iBeam,:),tempLatitude(iBeam,:),peakIonizationAlt,'linear','extrap');
            pfisrLongitude(iBeam,1) = interp1(tempAltitude(iBeam,:),tempLongitude(iBeam,:),peakIonizationAlt,'linear','extrap');  
        end
        plotData.latitude = repmat(pfisrLatitude,[1,size(plotData.zEnergyBin,2)]);
        plotData.longitude = repmat(pfisrLongitude,[1,size(plotData.zEnergyBin,2)]);
        plotData.projectionAltitude = peakIonizationAlt; % in KM
        plotData.timeNeutralAtmosphere = timeNeutralAtmosphere;
        
    case 'MagneticFieldMap'

        plotData.ionosphereCoord = readh5_variable_at_time(inputH5FileStr,...
            'ionosphereCoordGDZ',['/magneticMap/',p.Results.magFieldModelStr,'/'],[])';
        plotData.magFieldModelStr = p.Results.magFieldModelStr;
        plotData.time = unix_to_matlab_time(h5read(inputH5FileStr,['/magneticMap/',p.Results.magFieldModelStr,'/time']))';
    otherwise
        error('No or incorrect plotMode');
end



end
