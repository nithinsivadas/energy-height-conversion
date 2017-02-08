function [ B_abs, B_r, B_theta, B_eq, L, magneticLatitude, invariantLatitude ] =...
    get_dipole_mag_field( radiusInEarthRadii, polarAngle )
%% GET_DIPOLE_MAG_FIELD Dipole magnetic field parameters given input radius and polar angle
%-------------------------------------------------------------------------
% Input
%------
% radiusInEarthRadii: r/RE
% polarAngle        : theta (colatitude) [deg]
%-------------------------------------------------------------------------
% Output
%---------
% B_abs             : |B| of the dipole magnetic field [T]
% B_r               : Radial field of the dipole magnetic field [T]
% B_theta           : Field along polar angle of the dipole magnetic field [T]
% B_eq              : Total field at the equatorial plane
% L                 : The Mcllwain L-parameter L = radiusAtEquatorialPlane/RE of the
%                     particular field line
% magneticLatitude  : pi/2 - theta [deg]
% invariantLatitude : where a particular field line described by the L-shell
%                     intersects the ground
%----------------------------------------------------------------------------
% Modified: 25th Jan 2017 
% Created : 25th Jan 2017 
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%
r = radiusInEarthRadii;
theta = deg2rad(polarAngle);

% Mean magnetic field
B_mean = 3.12*10^-5; % [T] 

% Field along the radial direction
B_r = -2*B_mean*(r.^-3).*cos(theta);

% Field along theta
B_theta = -B_mean*(r.^-3).*sin(theta);

% Total field
B_abs = B_mean*(r.^-3).*(1+3*(cos(theta)).^2).^0.5;

% Converting theta (colatitude) into magnetic latitude
lambda = pi/2 - theta;
magneticLatitude = rad2deg(lambda);

% Estimating the L-shell we are on
L = r./(cos(lambda)).^2;
B_eq = B_mean/L.^3;

% Estimating the invariant latitude of that L-shell
Lambda = acos(L.^-0.5);
invariantLatitude = rad2deg(Lambda);

end

