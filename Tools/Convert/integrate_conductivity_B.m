function integratedData = integrate_conductivity_B(data,nBeams, projectionAltitude)
%% integrate_conductivitiy.m This function integrates conductivity
% along altitude
%--------------------------------------------------------------------------
% Input
%------
% data - sigma_P [(nhxnBeams) x nT]
%       i.e. ordered as beam1-alt1, beam1-alt2,..,beam1-altn, beam2-alt1 ... beam2-altn
%      - sigma_H [ "" ]
%      - time [nTx1]
%      - coords [(nhxnBeams) x 3] (lat,lon,alt)
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
sigma_H = zeros(nBeams,length(data.time));
sigma_P = zeros(nBeams,length(data.time));
lat = zeros(nBeams,1);
lon = zeros(nBeams,1);

for itime=1:1:length(data.time)
  for iBeam=1:1:nBeams
    altitudeGrid = data.conductanceCoordsB(iBeam,:,3);
    altitudeNo=find_altitude(altitudeGrid,projectionAltitude);
    lat(iBeam,1)=data.conductanceCoordsB(iBeam,altitudeNo,1);
    lon(iBeam,1)=data.conductanceCoordsB(iBeam,altitudeNo,2);
    
    sigma_H(iBeam,itime) = trapz(altitudeGrid,data.sigma_H_B(iBeam,:,itime))*1000; % km to m
    sigma_P(iBeam,itime) = trapz(altitudeGrid,data.sigma_P_B(iBeam,:,itime))*1000; % km to m
    
  end
end

integratedData.sigma_P = sigma_P;
integratedData.sigma_H = sigma_H;
integratedData.lat = lat;
integratedData.lon = lon;
integratedData.projectionAltitude = projectionAltitude;
integratedData.time = data.time;

end













%------ Note
% for itime=1:1:length(dataMag.time)
%     for iBeam=1:1:nBeams
%         avgLat      = mean(geodeticMagCoords(iBeam:nBeams:nData,1));
%         avgLon      = mean(geodeticMagCoords(iBeam:nBeams:nData,2));
%         altitude    = geodeticMagCoords(iBeam:nBeams:nData,3);
%         if iBeam<2
%             [altitudeNo] = find_altitude( altitude, 80);
%         end
%         altStart=iBeam+nBeams*(altitudeNo-1);
%
%         if itime==1
%             conductanceLat = geodeticMagCoords(altStart:nBeams:nData,1);
%             conductanceLon = geodeticMagCoords(altStart:nBeams:nData,2);
%             conductanceAlt = geodeticMagCoords(altStart:nBeams:nData,3);
%             nBeamConductance = length(conductanceAlt);
%             dataMag.conductanceCoords(1+(iBeam-1)*nBeamConductance:1:nBeamConductance*iBeam,1) = conductanceLat;
%             dataMag.conductanceCoords(1+(iBeam-1)*nBeamConductance:1:nBeamConductance*iBeam,2) = conductanceLon;
%             dataMag.conductanceCoords(1+(iBeam-1)*nBeamConductance:1:nBeamConductance*iBeam,3) = conductanceAlt;
%         end
%
%         conductanceData = get_conductivity_v1(conductanceAlt,...
%             dataMag.electronDensity(altStart:nBeams:nData,itime),avgLat,avgLon,dataMag.time(itime));
%
%         dataMag.sigma_H(1+(iBeam-1)*nBeamConductance:1:nBeamConductance*iBeam,itime) = conductanceData.sigma_H;
%         dataMag.sigma_P(1+(iBeam-1)*nBeamConductance:1:nBeamConductance*iBeam,itime) = conductanceData.sigma_P;
%
%
%     end
% end
