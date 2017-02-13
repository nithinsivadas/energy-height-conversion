function [differentialDirectionalNumFlux_latitude,...
    pitchAngleAtLatitude, invariantLatitude] =...
    get_pitch_angle_distribution...
    (differentialDirectionalNumFlux_equator,...
    pitchAngle, L, latitude )
%% get_pitch_angle_distribution.m Calculates the pitch angle distribution
% from input equatorial pitch angle distribution at any latitude on a 
% particular field line
%-------------------------------------------------------------------------
% Input
%------
% differenitalDirectionalNumFlux_equator : j(alpha, E) [1 / (m sr s eV)] 
% pitchAngle                             : [deg]
% L                                      : magnetic L-shell (r_eq/RE)
% latitude                               : [deg]
%-------------------------------------------------------------------------
% Output
%-------
% differentialDirectionalNumFlux_latitude : j(alpha_lambda, E) [1/ m sr s eV)]
% pitchAngleAtLatitude                    : Pitch angle associated with
%                                           values of diff num flux
% invariantLatitude                       : The invariant latitude of the
%                                           point
%%

% Initializing the Equatorial Pitch Angle Distribution
j_eq = differentialDirectionalNumFlux_equator;
alpha_eq = deg2rad(pitchAngle);

% Calculating the magnetic field dipole parameters
lambda = deg2rad(latitude);
theta = pi/2-lambda;
radiusInEarthRadii = L*(cos(lambda)).^2;

[B_lambda,B_r,B_theta,B_eq,L,magneticLatitude,invariantLatitude]=...
    get_dipole_mag_field( radiusInEarthRadii, rad2deg(theta) );

% Estimating the pitch angle at the non-equatorial point along the magnetic
% field line
alpha_lambda = deg2rad(0:1:180);
alpha_eq_new = asin(sin(alpha_lambda).*(B_eq./B_lambda).^0.5);

j_lambda = interp1(alpha_eq, j_eq, alpha_eq_new,'linear','extrap');

pitchAngleAtLatitude = rad2deg(alpha_lambda);
differentialDirectionalNumFlux_latitude = j_lambda;

%%
%----------------------------------------------------------------------------
% Modified: 8th Feb 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
end

