function integratedData = integrate_conductivity(data,nBeams, projectionAltitude)
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

ndata = (length(data.conductanceCoords(:,1)))/nBeams;

for itime=1:1:length(data.time)
  for iBeam=1:1:nBeams
    sigma_H(iBeam,itime) = 0;
    sigma_P(iBeam,itime) = 0;
    altitudeGrid = data.conductanceCoords(1+(iBeam-1)*ndata:1:iBeam*ndata,3);
    altitudeNo=find_altitude(altitudeGrid,projectionAltitude);
    lat(iBeam,1)=data.conductanceCoords(altitudeNo + (iBeam-1)*ndata,1);
    lon(iBeam,1)=data.conductanceCoords(altitudeNo + (iBeam-1)*ndata,2);
    for idata=1:1:ndata-1
      deltaAltitude = (data.conductanceCoords(idata+1 +(iBeam-1)*ndata,3)...
      - data.conductanceCoords(idata +(iBeam-1)*ndata,3))*1000;  % converting km to m
      sigma_H(iBeam,itime) = sigma_H(iBeam,itime)...
      + data.sigma_H(idata + (iBeam-1)*ndata,itime)*deltaAltitude;
      sigma_P(iBeam,itime) = sigma_P(iBeam,itime)...
      + data.sigma_P(idata + (iBeam-1)*ndata,itime)*deltaAltitude;
    end
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
