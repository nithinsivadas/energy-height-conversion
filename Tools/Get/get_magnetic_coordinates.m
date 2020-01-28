function [MLAT, MLON, MLT] = get_magnetic_coordinates(GDZ,time)
%get_magnetic_coordinates Calculate magnetic coordinates (MLAT, MLON, MLT) given a
%particular geodetic coordinate and time. 
% Input
%     GDZ  - [Alt, Lat, Lon] in [km, deg, deg]
% Output
%     MLAT - Magnetic latitude in degrees
%     MLON - Magnetic longitude in degrees
%     MLT  - Magnetic Local Time in hours
%--------------------------------------------------------------
% Created on: 30th May 2019
% Updated on: 11th Jan 2020; Error in the previous code, MAG not SM. 
% Created by: Nithin Sivadas
% Comments:
% Requires IRBEM Coordinate transform
%--------------------------------------------------------------

% xSM = onera_desp_lib_coord_trans(GDZ,[0,4],time(:));
% xSM_sph = onera_desp_lib_rotate(xSM,'car2sph');
% % [az, el, r] = cart2sph(xSM(1),xSM(2), xSM(3));
% MLON = convert_longitude(xSM_sph(:,3)+90,'360to180');
% MLT = (xSM_sph(:,3)+180)*(24/360);
% MLAT = xSM_sph(:,2);

% Laundal & Richmond 2017
% MLT2 = (time(:)-floor(time(:))).*24 + (xSM_sph(:,3)+90 -72.63)./15;

xMAG = onera_desp_lib_coord_trans(GDZ,[0,6],time(:));
MAG_sph = onera_desp_lib_rotate(xMAG,'car2sph');
MLON = convert_longitude(MAG_sph(:,3),'360to180');
MLT = mod((time-floor(time)).*24 + mod(MAG_sph(:,3)-72.63,360)./15,24);
MLAT = MAG_sph(:,2);

end

