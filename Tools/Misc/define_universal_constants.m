function [C] = define_universal_constants()
%define_universal_constants.m Defining universal constants
%   Output
%   C.info  - contains the name of the constant and units
%   C.c     - speed of light`
%   C.kb    - boltazmann constant
%   C.e     - elementary charge
%   C.epsilon0 - electric constant

  C.c = 299792458; %m/s
  C.info.c = 'Speed of light [m/s]';

  C.kb = 1.38064852e-23;
  C.info.kb = 'Boltzmann Constant [m^2 kg s^-2 K^-1]';

  C.e = 1.60217662e-19;
  C.info.e = 'Elementary charge [coulombs]';

  C.epsilon0 = 8.85418782e-12;
  C.info.epsilon0 = 'Electric constant [m-3 kg-1 s4 A2]';

   C.RE = 6.371*10^6;
   C.info.RE = 'Average Radius of Earth [m]';
end
