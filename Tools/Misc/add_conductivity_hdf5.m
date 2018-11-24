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

% Max time instances that can be recorded to HDF5 file at once
% to avoid memory errors
nTimeMax = 1000;

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
nTime = length(time);

% nTime = 3;
% pTimeIndx = find_time(time,'26-Mar-2008 11:08');

% Calculating conductivity
multiWaitbar('Calculating conductivity...',0);
di = 1./nTime;
for iTime = 1:1:nTime
    alt1 = linspace(min(alt(:)),max(alt(:)),2*length(alt(:,1)));
    lat1 = lat*ones(size(alt1));
    lon1 = lon*ones(size(alt1));
    msisData = msis_irbem(time(iTime), [alt1',lat1',lon1']);
   for iBeam = 1:1:nBeams
       
        inputNe = squeeze(Ne(:,iBeam,iTime));
        if setInterpNe
        inputNe(inputNe<0) = nan;
        inputNe = interp_nans(inputNe);
        end
       inputAlt = squeeze(alt(:,iBeam));
       
       % Compromise to increase the code-speed by 10 times. 
       inputMsisData.AltTemp = interp1(alt1,msisData.AltTemp,inputAlt);
       inputMsisData.O = interp1(alt1,msisData.O,inputAlt);
       inputMsisData.N2 = interp1(alt1,msisData.N2,inputAlt);
       inputMsisData.O2 = interp1(alt1,msisData.O2,inputAlt);
       inputMsisData.totalNumberDensity = interp1(alt1,msisData.totalNumberDensity,inputAlt);
       
       [data, input] = get_conductivity_v2(inputAlt, inputNe, lat, lon,...
           time(iTime), 0, {'all'}, false,[],inputMsisData);
       
       out.sigma_P(:,iBeam,iTime) = real(data.pedersonConductivity);
       out.sigma_H(:,iBeam,iTime) = real(data.hallConductivity);
       out.alt(:,iBeam,1) = data.altitude;
       out.time(iTime) = posixtime(datetime(time(iTime),'ConvertFrom','datenum'));
       in.Ne(:,iBeam,iTime) = input.electronDensity;
%        in.Ne0(:,iBeam,iTime) = inputNe;
   end
   multiWaitbar('Calculating conductivity...','Increment',di);
end

multiWaitbar('Writing to HDF5...',0);

for iTimeMax = 1:1:ceil(nTime./nTimeMax)
    multiWaitbar('Writing to HDF5...',0.5);
    timeIndx = 1 + (iTimeMax-1)*nTimeMax:1:min(iTimeMax*nTimeMax,nTime);
    write_h5_dataset(outputH5Str,'/conductivity/time',...
        time(timeIndx)',1);
    write_h5_dataset(outputH5Str,'/conductivity/sigmaP',...
        permute(out.sigma_P(:,:,timeIndx),[3 2 1]),1);
    write_h5_dataset(outputH5Str,'/conductivity/sigmaH',...
        permute(out.sigma_H(:,:,timeIndx),[3 2 1]),1);
    write_h5_dataset(outputH5Str,'/conductivity/inputNe',...
        permute(in.Ne(:,:,timeIndx),[3 2 1]),1);
end
write_h5_dataset(outputH5Str,'/conductivity/alt',...
        out.alt',0);

write_h5_dataset_attribute(outputH5Str,'/conductivity/time',...
    [],[],[],'time');
write_h5_dataset_attribute(outputH5Str,'/conductivity/sigmaP',...
    'Peserson conductivity profiles','[nTime x nBeams x nAlt]','[S m^-^1]');
write_h5_dataset_attribute(outputH5Str,'/conductivity/sigmaH',...
    'Hall conductivity profiles','[nTime x nBeams x nAlt]','[S m^-^1]');
write_h5_dataset_attribute(outputH5Str,'/conductivity/inputNe',...
    'Electron density profile used to calculate conductivity','[nTime x nBeams x nAlt]','[m^-^3]');
write_h5_dataset_attribute(outputH5Str,'/conductivity/alt',...
    'Altitude','[nBeams x nAlt]','[km]');
multiWaitbar('Writing to HDF5...',1);

end
