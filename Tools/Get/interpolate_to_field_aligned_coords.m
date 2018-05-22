function [amisrData] = interpolate_to_field_aligned_coords(amisrData,timeMinStr,timeMaxStr)
%interpolate_to_field_aligned_coords Summary of this function goes here
% Input
%   amisrData.origCartCoords
%   amisrData.magCartCoords
% Input
%   amisrData.magElectronDensity

% Error check
if ~isfield(amisrData,'origCartCoords')
    error(['Original coordinates field (origCartCoords)',...
        'not available in the structure amisrData']);
end
if ~isfield(amisrData,'magCartCoords')
    error(['Magnetic field aligned coordinates field (magCartCoords)',...
        'not available in the structure amisrData']);
end 

if nargin < 3 || isempty(timeMaxStr)
    timeMaxStr = datestr(max(amisrData.time(1,:)));
end
if nargin < 2 || isempty(timeMinStr)
    timeMinStr = datestr(min(amisrData.time(1,:)));
end

timeMinIndx = find_time(amisrData.time(1,:),timeMinStr);
timeMaxIndx = find_time(amisrData.time(1,:),timeMaxStr);
dtime = 1./(length(timeMinIndx:1:timeMaxIndx));
multiWaitbar('Interpolate Ne along field-aligned coords',0);
    for itime = timeMinIndx:timeMaxIndx
        multiWaitbar('Interpolate Ne along field-aligned coords','Increment',dtime);
        ne=amisrData.electronDensity(:,:,itime);
        F = scatteredInterpolant(amisrData.origCartCoords.xEast(:),...
            amisrData.origCartCoords.yNorth(:),...
            amisrData.origCartCoords.zUp(:),ne(:),'linear','none');
        amisrData.magElectronDensity(:,:,itime) = ...
            F(amisrData.magCartCoords.xEast,amisrData.magCartCoords.yNorth,...
            amisrData.magCartCoords.zUp);
        
%       Removing the NANs of Ne
        dne = amisrData.dNeFrac(:,:,itime);
        dneModel = (nanmean(amisrData.dNeFrac,3));
        dne(isnan(dne))=dneModel(isnan(dne));
        
        dF = scatteredInterpolant(amisrData.origCartCoords.xEast(:),...
            amisrData.origCartCoords.yNorth(:),...
            amisrData.origCartCoords.zUp(:),dne(:),'linear','none');
        amisrData.magdNeFrac(:,:,itime) = ...
            dF(amisrData.magCartCoords.xEast,amisrData.magCartCoords.yNorth,...
            amisrData.magCartCoords.zUp);
    end
% multiWaitbar('Interpolate Ne along field-aligned coords','Close');
end

