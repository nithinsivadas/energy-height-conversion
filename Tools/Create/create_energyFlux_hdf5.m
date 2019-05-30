function [data] = create_energyFlux_hdf5(dataStructure,amisrMagDataStructure,h5FileStr,mode)
%create_energyFlux_hdf5 Creates/adds energy flux from inverted PFISR measurements
% of electron density profiles to hdf5 files. The energy flux is input in the form
% of a data structure (which is output of another funciton).
%   Detailed explanation goes here
%   if Mode = 'dataInv': then the structure is in the following format
%              dataInv(nBeam).energyFLux(nEnergy/nAlt,nTime);
if nargin<4
    mode = 'dataInv';
end
if nargin<3
    h5FileStr = 'output.h5';
end

if strcmp(mode,'dataInv')
    %% Formatting data
    nBeams=length(dataStructure);
    fieldNames = fieldnames(dataStructure(1));
    nFields = length(fieldNames);
    nTime = length(dataStructure(1).time);
    nAlt = length(dataStructure(1).alt);
    nEnergy = length(dataStructure(1).energyBin);

    for iField = 1:1:nFields
        data{1,iField} = fieldNames{iField};
        if iField == 1|| iField == 2 || iField == 3 || iField == 10
            for iBeam = 1:1:nBeams
            data{2,iField}(:,iBeam,:)=dataStructure(iBeam).(fieldNames{iField})';
            end
        elseif iField == 12 || iField == 11 || iField == 5
            for iBeam = 1:1:nBeams
            data{2,iField}(:,iBeam,:)=dataStructure(iBeam).(fieldNames{iField});
            end
        elseif iField == 4 || iField == 6
            for iBeam = 1:1:nBeams
            data{2,iField}(:,iBeam)=dataStructure(iBeam).(fieldNames{iField});
            end
        else
            data{2,iField}=dataStructure(1).(fieldNames{iField})';
        end
    end

    inputFieldNames = fieldnames(amisrMagDataStructure);
    nInputFields = length(inputFieldNames);

    for iField=1:1:nInputFields
        data{1,iField+nFields} = inputFieldNames{iField};
            if iField+nFields==17
                data{2,iField+nFields}=...
                        amisrMagDataStructure(1).(inputFieldNames{iField});
            else
                for iBeam = 1:1:nBeams
                    data{2,iField+nFields}(:,iBeam,:)=...
                        amisrMagDataStructure(iBeam).(inputFieldNames{iField});
                end
            end
    end

    data{3,1} = 'nTime x nBeams x nEnergyBin';
    data{4,1} = 'Differential number flux';
    data{5,1} = '[m^-2 s^-1 eV^-1]';

    data{3,2} = 'nTime x nBeams x nEnergyBin';
    data{4,2} = 'Differential energy flux';
    data{5,2} = '[eV m^-2 sr^-1 s^-1 eV^-1]';

    data{3,3} = 'nTime x nBeams x nEnergyBin';
    data{4,3} = 'Error in differential energy flux';
    data{5,3} = '[eV m^-2 sr^-1 s^-1 eV^-1]';

    data{3,4} = 'nTime x nBeams';
    data{4,4} = 'Reduced chi-squared value of MEM Inversion Fit sqrt(sum((data-fit)/error)^2)/DOF';
    data{5,4} = '[a.u.]';

    data{3,5} = 'nTime x nBeams x nAltitude';
    data{4,5} = 'Production rate derived from the estimated flux - A*flux ';
    data{5,5} = '[m^-3 s^-1]';

    data{3,6} = 'nTime x nBeams';
    data{4,6} = 'The maximum number of iterations before convergence';
    data{5,6} = '[number]';

    data{3,7} = '1 x nEnergyBins';
    data{4,7} = 'The energy values of each energy-spectral bin';
    data{5,7} = '[eV]';

    data{3,8} = 'nTime x 1';
    data{4,8} = 'DateTime Vector';
    data{5,8} = '[matlab units]';

    data{3,9} = '1 x nAltitude';
    data{4,9} = 'Altitude of Ne measurement';
    data{5,9} = '[km]';

    data{3,10} = 'nEnergy x nBeams x nHeight';
    data{4,10} = 'Production rate per unit differential number flux';
    data{5,10} = '[m^-1 eV^-1]';

    data{3,11} = 'nTime x nBeams x nAltitude';
    data{4,11} = 'The production rate which was calculated from Ne and D-region chemistry model Vickrey et al., 1982';
    data{5,11} = '[m^-3 s^-1]';

    data{3,12} = 'nTime x nBeams x nAltitude';
    data{4,12} = 'The error in production rate';
    data{5,12} = '[m^-3 s^-1]';

    data{3,13} = 'nTime x nBeams x nAltitude';
    data{4,13} = 'Magnetically field aligned electron density';
    data{5,13} = '[m^-3]';

    data{3,14} = 'nTime x nBeams x nAltitude';
    data{4,14} = 'Fractional error in electron density';
    data{5,14} = '[ratio]';

    data{3,15} = 'nCoords x nBeams x nAltitude';
    data{4,15} = 'Magnetically field aligned Geodetic coordinates [Lat, Lon, Alt]';
    data{5,15} = '[deg, deg, km]';

    data{3,16} = 'nCoords x nBeams x nAltitude';
    data{4,16} = 'Magnetically field aligned Cartesian coordinates [xEast, yNorth, Alt]';
    data{5,16} = '[km, km, km]';

    data{3,17} = '1x1';
    data{4,17} = 'Projection Altitude; the altitude from which the mag. field aligned coordinates are drawn';
    data{5,17} = '[km]';

    %% Writing into HDF5 File
    dFields = 1./(nFields+nInputFields);
    multiWaitbar('Write energyflux into HDF5',0);
    for iField=1:1:nFields
        multiWaitbar('Write energyflux into HDF5','Increment',dFields);
        if iField ~= 6 && iField ~= 7 && iField ~= 8 && iField ~= 9 && iField ~= 4
            h5create(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],size(permute(data{2,iField},[3 2 1])));
            h5write(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],permute(data{2,iField},[3 2 1]));
        else
            h5create(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],size((data{2,iField})'));
            h5write(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],(data{2,iField})');
        end
        h5writeatt(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],'Dimensions',data{3,iField});
        h5writeatt(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],'Description',data{4,iField});
        h5writeatt(h5FileStr,['/energyFluxFromMaxEnt/',data{1,iField}],'Units',data{5,iField});
    end

    for iField=nFields+1:nFields+nInputFields
        multiWaitbar('Write energyflux into HDF5','Increment',dFields);
        if iField == 13 || iField == 14
            h5create(h5FileStr,['/inputData/',data{1,iField}],size(permute(data{2,iField},[3 2 1])));
            h5write(h5FileStr,['/inputData/',data{1,iField}],permute(data{2,iField},[3 2 1]));
            h5writeatt(h5FileStr,['/inputData/',data{1,iField}],'Dimensions',data{3,iField});
            h5writeatt(h5FileStr,['/inputData/',data{1,iField}],'Description',data{4,iField});
            h5writeatt(h5FileStr,['/inputData/',data{1,iField}],'Units',data{5,iField});
        elseif iField == 17
            h5create(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],size((data{2,iField})));
            h5write(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],(data{2,iField}));
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Dimensions',data{3,iField});
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Description',data{4,iField});
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Units',data{5,iField});

        else
            h5create(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],size(permute(data{2,iField},[3 2 1])));
            h5write(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],permute(data{2,iField},[3 2 1]));
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Dimensions',data{3,iField});
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Description',data{4,iField});
            h5writeatt(h5FileStr,['/magneticFieldAlignedCoordinates/',data{1,iField}],'Units',data{5,iField});
        end
    end
    h5writeatt(h5FileStr,'/','creation_date',datestr(now));
    h5writeatt(h5FileStr,'/','duration_contained_in_file',[datestr(data{2,8}(1)),' - ',datestr(data{2,8}(end))]);
%     multiWaitbar('Write energyflux into HDF5','Close');
end

end
