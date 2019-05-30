function plotData = get_2D_plot_inputs_time_independent(inputH5FileStr,varargin)
% get_2D_plot_inputs_at_time.m Get 'OpticalImage', 'EnergyFluxMap',
% 'MagneticFieldMap' data form hdf5 file that is time independent
%   Detailed explanation goes here

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

expectedMaps = {'OpticalImage','EnergyFluxMap','MagneticFieldMap','NoMap'};
expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};
expectedSites = {'gako','fykn','mcgr','whit','inuv','kian','dasc','pokerFlat'};

addParameter(p,'plotModeStr','NoMap',@(x) any(validatestring(x,expectedMaps)));
addParameter(p,'site','pokerFlat',@(x) any(validatestring(x,expectedSites)));
addParameter(p,'magFieldModelStr','TS96',@(x) any(strcmp(x,expectedMagFieldModels)));
addParameter(p,'timeNeutralAtmosphere',nan,validScalarPosNum);
addParameter(p,'energySlice',100,validScalarPosNum); %In keV
addParameter(p,'peakIonizationAltitude',nan); %In Km
addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));

parse(p,inputH5FileStr,varargin{:});


if nargin < 1
    error('inputH5FileStr not specified');
end


switch p.Results.plotModeStr
    case 'OpticalImage'
        if strcmp(p.Results.site,'pokerFlat')
            siteStr = 'dasc';
        else
            siteStr = p.Results.site;
        end
        siteStr = ['/',upper(siteStr),'/'];
        plotData.latitude = readh5_variable_at_time(inputH5FileStr,'lat',...
            siteStr,[]);
        plotData.longitude = readh5_variable_at_time(inputH5FileStr,'lon',...
            siteStr,[]);
        plotData.azimuth = readh5_variable_at_time(inputH5FileStr,'az',...
            siteStr,[]);
        plotData.elevation = readh5_variable_at_time(inputH5FileStr,'el',...
            siteStr,[]);
        plotData.time = unix_to_matlab_time(h5read(inputH5FileStr,['/',upper(siteStr),'/time']))';
        plotData.background = h5read(inputH5FileStr,['/',upper(siteStr),'/background']);
        plotData.sensorloc = h5read(inputH5FileStr,['/',upper(siteStr),'/sensorloc']);
        if ~strcmp(siteStr,'/DASC/')
            plotData.sensorloc(2) = convert_longitude(plotData.sensorloc(2),'360to180');
            plotData.longitude = convert_longitude(plotData.longitude,'360to180');
        end
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

        if isnan(p.Results.peakIonizationAltitude)

            if isnan(p.Results.timeNeutralAtmosphere)
                timeNeutralAtmosphere=median(plotData.time);
            else
                timeNeutralAtmosphere=p.Resilts.timeNeutralAtmosphere;
            end

            for iEnergy = 1:1:size(plotData.zEnergyBin,2)
            peakIonizationAlt(iEnergy) = calculate_peak_altitude_of_ionization(plotData.zEnergyBin(1,iEnergy),...
                timeNeutralAtmosphere,tempLatitude(1,:)',tempLongitude(1,:)',...
                tempAltitude(1,:)');
            end
            fprintf(['PFISR Energy Flux Map projected at ',num2str(peakIonizationAlt),' km\n']);
            nBeams = size(tempLatitude,1);
            for iEnergy = 1:1:size(plotData.zEnergyBin,2)
                for iBeam = 1:1:nBeams
                    pfisrLatitude(iBeam,iEnergy) = interp1(tempAltitude(iBeam,:),tempLatitude(iBeam,:),peakIonizationAlt(iEnergy),'linear','extrap');
                    pfisrLongitude(iBeam,iEnergy) = interp1(tempAltitude(iBeam,:),tempLongitude(iBeam,:),peakIonizationAlt(iEnergy),'linear','extrap');
                end
            end
        else
            timeNeutralAtmosphere = nan;
            peakIonizationAlt = p.Results.peakIonizationAltitude*ones(1,size(plotData.zEnergyBin,2));

            nBeams = size(tempLatitude,1);
            iEnergy = 1;
            for iBeam = 1:1:nBeams
                pfisrLatitude(iBeam,iEnergy) = interp1(tempAltitude(iBeam,:),tempLatitude(iBeam,:),peakIonizationAlt(iEnergy),'linear','extrap');
                pfisrLongitude(iBeam,iEnergy) = interp1(tempAltitude(iBeam,:),tempLongitude(iBeam,:),peakIonizationAlt(iEnergy),'linear','extrap');
            end
            pfisrLatitude = repmat(pfisrLatitude,1,size(plotData.zEnergyBin,2));
            pfisrLongitude = repmat(pfisrLongitude,1,size(plotData.zEnergyBin,2));

        end

%         plotData.latitude = repmat(pfisrLatitude,[1,size(plotData.zEnergyBin,2)]);
%         plotData.longitude = repmat(pfisrLongitude,[1,size(plotData.zEnergyBin,2)]);
        plotData.latitude = pfisrLatitude;
        plotData.longitude = pfisrLongitude;
        plotData.projectionAltitude = peakIonizationAlt; % in KM
        plotData.timeNeutralAtmosphere = timeNeutralAtmosphere;

    case 'MagneticFieldMap'
        plotData.ionosphereCoord = readh5_variable_at_time(inputH5FileStr,...
            'ionosphereCoordGDZ',['/magneticMap/',p.Results.magFieldModelStr,'/'],[])';
        plotData.time = unix_to_matlab_time(h5read(inputH5FileStr,['/magneticMap/',p.Results.magFieldModelStr,'/time']))';
    otherwise
        error('No or incorrect plotMode');
end



end
