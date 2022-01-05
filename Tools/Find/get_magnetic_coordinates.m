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
% Created by: Nithin Sivadas
% Comments:
% Requires IRBEM Coordinate transform
%--------------------------------------------------------------

xSM = onera_desp_lib_coord_trans(GDZ,[0,4],time);
xSM_sph = onera_desp_lib_rotate(xSM,'car2sph');
% [az, el, r] = cart2sph(xSM(1),xSM(2), xSM(3));
MLON = convert_longitude(xSM_sph(3)+90,'360to180');
MLT = (xSM_sph(3)+180)*(24/360);
MLAT = xSM_sph(2);


end

