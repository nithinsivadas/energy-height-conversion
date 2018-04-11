function [msisData,F107A,F107,APH] = msis_irbem(timeMSIS,coords,msisVersionStr,coordType)
%UNTITLED2 Summary of this function goes here
%% Coords - [alt; lat; lon]
%   Detailed explanation goes here
if nargin<4
    coordType=0;
end
if nargin<3
    msisVersionStr='nrlmsise00';
end

timeaph = min(timeMSIS(:));
[F107A, F107, APH] = f107_aph(timeaph);
coords=reshape(coords,3,[]);
alt = coords(1,:);
lat = coords(2,:);
lon = coords(3,:);
nCoords=length(alt);
msisData = onera_desp_lib_msis(msisVersionStr,timeMSIS,[alt(:),lat(:),lon(:)],...
    coordType,repmat(F107A,nCoords,1),repmat(F107,nCoords,1),repmat(APH,nCoords,1));
end

