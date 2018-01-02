function integratedData = integrate_conductivity(data,nBeams, projectionAltitude, altRange)
%% integrate_conductivitiy.m This function integrates conductivity
% along altitude, along magneticfield aligned beams
%--------------------------------------------------------------------------
% Input
%------
% data      - sigma_P_B [(nBeams x nh) x nT]
%             calculated along mangetic field aligned nBeams
%           - sigma_H_B [(nBeams x nh) x nT]
%             calculated along magneticfield aligned beams
%           - time [nTx1]
%           - conductanceCoordsB [(nBeams x nh) x 3] (lat,lon,alt)
%             calculated along magneticfield aligned beams
% nBeams    - total number of radar beams
% projectionAltitude - Altitude [km] at which the map ought to be projected
% altRange  - [minAlt maxAlt] the range of altitude to be used for integration

%--------------------------------------------------------------------------
% Output
%-------
% integratedData - sigma_H [nBeams x nT]
%                  sigma_P [nBeams x nT]
%                  lat     [nBeams x 1]
%                  lon     [nBeams x 1]
%                  projectionAltitude [1]
%                  time   [nT x 1]
%                  ratio  [nBeams x nT] : ratio of Hall to Pedersen conductance

%--------------------------------------------------------------------------
% Modified: 2nd Jan 2018
% Created : 31st Nov 2017
% Author  : Nithin Sivadas
% Ref     :
%--------

  if nargin<4
    altRange = [min(data.conductanceCoordsB(1,:,3)) max(data.conductanceCoordsB(1,:,3))];
  end

 conductanceCoordsA = data.conductanceCoordsB;
 sigma_H_A = data.sigma_H_B;
 sigma_P_A = data.sigma_P_B;

  % Initializing conductance & position variables
  sigma_H = zeros(nBeams,length(data.time));
  sigma_P = zeros(nBeams,length(data.time));
  lat = zeros(nBeams,1);
  lon = zeros(nBeams,1);

  for itime=1:1:length(data.time)
    for iBeam=1:1:nBeams
      % Recording appropriate position beased on projectionAltitude
      minAltNo=find_altitude(conductanceCoordsA(iBeam,:,3),altRange(1));
      maxAltNo=find_altitude(conductanceCoordsA(iBeam,:,3),altRange(2));
      altitudeGrid = conductanceCoordsA(iBeam,minAltNo:maxAltNo,3);
      altitudeNo=find_altitude(altitudeGrid,projectionAltitude);
      lat(iBeam,1)=conductanceCoordsA(iBeam,altitudeNo,1);
      lon(iBeam,1)=conductanceCoordsA(iBeam,altitudeNo,2);
      % Integrating conductivity along each beam
      sigma_H(iBeam,itime) = trapz(altitudeGrid,sigma_H_A(iBeam,minAltNo:maxAltNo,itime))*1000; % km to m
      sigma_P(iBeam,itime) = trapz(altitudeGrid,sigma_P_A(iBeam,minAltNo:maxAltNo,itime))*1000; % km to m
      % Is it okay to integrate along altitude, and not along range?
    end
  end

  % Storing output variables
  integratedData.sigma_P = sigma_P;
  integratedData.sigma_H = sigma_H;
  integratedData.lat = lat;
  integratedData.lon = lon;
  integratedData.projectionAltitude = projectionAltitude;
  integratedData.time = data.time;
  integratedData.ratio = sigma_H./sigma_P;
end


% A different version of the function with sigma_P [ndata x nT]
%function integratedData = integrate_conductivity(data,nBeams, projectionAltitude)
%% integrate_conductivitiy.m This function integrates conductivity
% along altitude
%--------------------------------------------------------------------------
% Input
%------
% data - sigma_P [(nhxnBeams) x nT]
%       i.e. ordered as beam1-alt1, beam1-alt2,..,beam1-altn, beam2-alt1 ... beam2-altn
%      - sigma_H [ "" ]
%      - time [nTx1]
%      - conductanceCoords [(nhxnBeams) x 3] (lat,lon,alt)
%--------------------------------------------------------------------------
% Output
%-------
% integratedData - sigma_H [nbeams x nT]
%--------------------------------------------------------------------------
% Modified: 25th Sep 2016
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%--------

% ndata = (length(data.conductanceCoords(:,1)))/nBeams;
%
% for itime=1:1:length(data.time)
%   for iBeam=1:1:nBeams
%     sigma_H(iBeam,itime) = 0;
%     sigma_P(iBeam,itime) = 0;
%     altitudeGrid = data.conductanceCoords(1+(iBeam-1)*ndata:1:iBeam*ndata,3);
%     altitudeNo=find_altitude(altitudeGrid,projectionAltitude);
%     lat(iBeam,1)=data.conductanceCoords(altitudeNo + (iBeam-1)*ndata,1);
%     lon(iBeam,1)=data.conductanceCoords(altitudeNo + (iBeam-1)*ndata,2);
%     for idata=1:1:ndata-1
%       deltaAltitude = (data.conductanceCoords(idata+1 +(iBeam-1)*ndata,3)...
%       - data.conductanceCoords(idata +(iBeam-1)*ndata,3))*1000;  % converting km to m
%       sigma_H(iBeam,itime) = sigma_H(iBeam,itime)...
%       + data.sigma_H(idata + (iBeam-1)*ndata,itime)*deltaAltitude;
%       sigma_P(iBeam,itime) = sigma_P(iBeam,itime)...
%       + data.sigma_P(idata + (iBeam-1)*ndata,itime)*deltaAltitude;
%     end
%   end
% end
%
% integratedData.sigma_P = sigma_P;
% integratedData.sigma_H = sigma_H;
% integratedData.lat = lat;
% integratedData.lon = lon;
% integratedData.projectionAltitude = projectionAltitude;
% integratedData.time = data.time;
%
% end
