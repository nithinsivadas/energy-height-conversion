function [differentialDirectionalNumFlux_latitude,...
    pitchAngleAtLatitude, invariantLatitude] =...
    get_pitch_angle_distribution...
    (differentialDirectionalNumFlux_equator,...
    pitchAngle, L, latitude )
% get_pitch_angle_distribution Calculates the pitch angle distribution
% from input equatorial pitch angle distribution at any latitude on a 
% particular field line

% Input
% differenitalDirectionalNumFlux_equator : j(alpha, E) [1 / (m sr s eV)] 
% pitchAngle : [deg]
% L          : magnetic L-shell (r_eq/RE)
% latitude   : [deg]

j_eq = differentialDirectionalNumFlux_equator;
alpha_eq = deg2rad(pitchAngle);

lambda = deg2rad(latitude);
theta = pi/2-lambda;
radiusInEarthRadii = L*(cos(lambda)).^2;

[B_lambda,B_r,B_theta,B_eq,L,magneticLatitude,invariantLatitude]=...
    get_dipole_mag_field( radiusInEarthRadii, rad2deg(theta) );

alpha_lambda = deg2rad(0:1:180);
alpha_eq_new = asin(sin(alpha_lambda).*(B_eq./B_lambda).^0.5);

j_lambda = interp1(alpha_eq, j_eq, alpha_eq_new,'linear','extrap');

pitchAngleAtLatitude = rad2deg(alpha_lambda);
differentialDirectionalNumFlux_latitude = j_lambda;

end

