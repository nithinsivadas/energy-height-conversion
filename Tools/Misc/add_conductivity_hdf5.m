function [out, in] = add_conductivity_hdf5(inputH5Str,outputH5Str,setInterpNe)
%add_conductivity_hdf5.m calculates Pedersen and Hall conductivities, from
%magnetically field aligned electron densities obtained from the
%inputH5Str, and stores it in outputH5Str. 
%-------------------------------------------------------------------------
% Input
%------
%       inputH5Str - Filename of the input HDF5 file, that contains
%                   /inputData/Ne
%                   /energyFluxFromMaxEnt/time
%                   /magneticFieldAlignedCoordinates/magGeodeticLatLonAlt/
%-------------------------------------------------------------------------
% Output
%---------
%       
%   outputArguments.sigma_P  - Pederson conductivity [S/m]
%   outputArguments.sigma_H  - Hall conductivity [S/m]
%   outputArguments.alt      - Altitude [km]
%   and  a lot of other things
%----------------------------------------------------------------------------
% Modified: 20th Nov 2018,12th Jan 2018
% Created : 27th Sep 2017
% Author  : Nithin Sivadas
% Ref     :
% Comments: 
% 
%----------------------------------------------------------------------------

% Extracting Data from InputH5File
time = h5read(inputH5Str,'/energyFluxFromMaxEnt/time'); 
Ne = h5read(inputH5Str,'/inputData/Ne');
coords = h5read(inputH5Str,...
    '/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt');

% Converting Matrices to arrays, to avoid forloop
alt = squeeze(coords(:,:,3));
lat = mean(reshape(squeeze(coords(:,:,1)),1,[]));
lon = mean(reshape(squeeze(coords(:,:,2)),1,[]));

% Numbers necessary to convert the array back to matrix
nBeams = size(coords,2);
nAlt = size(coords,1);
nTime = length(time);
pTimeIndx = find_time(time,'26-Mar-2008 11:08');
% Calculating conductivity
tic
for iTime = pTimeIndx:1:pTimeIndx
   for iBeam = 1:1:nBeams
       
        inputNe = squeeze(Ne(:,iBeam,iTime));
        if setInterpNe
        inputNe(inputNe<0) = nan;
        inputNe = interp_nans(inputNe);
        end
       inputAlt = squeeze(alt(:,iBeam));
       
       [data, input] = get_conductivity_v2(inputAlt, inputNe, lat, lon,...
           time(iTime), 0, {'all'}, false);
       
       out.sigma_P(:,iBeam,iTime) = data.pedersonConductivity;
       out.sigma_H(:,iBeam,iTime) = data.hallConductivity;
       out.alt(:,iBeam,iTime) = data.altitude;
       out.time(iTime) = time(iTime);
       in.Ne(:,iBeam,iTime) = input.electronDensity;
       in.Ne0(:,iBeam,iTime) = inputNe;
   end
end
toc
figure; 
semilogx(squeeze(out.sigma_P(:,13,iTime)),squeeze(out.alt(:,13,iTime))); 
hold on; 
semilogx(squeeze(out.sigma_H(:,13,iTime)),squeeze(out.alt(:,13,iTime)),'-');
title(datestr(time(iTime)));
% alt1 = convert_array_to_3D(alt,nAlt,nBeams); %Converting

end

function matrix = convert_array_to_3D(array,nAlt,nBeams)
matrix = reshape(array,nAlt,nBeams); %Converting
end

function array = convert_matrix_to_1D(matrix)
array = reshape(matrix,1,[]); %Converting
end

