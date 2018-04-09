function [C] = define_universal_constants()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
C.c = 299792458; %m/s
C.info.c = 'Speed of light [m/s]';
C.kb = 1.38064852e-23;
C.info.kb = 'Boltzmann Constant [m^2 kg s^-2 K^-1]';
C.e = 1.60217662e-19;
C.info.e = 'Elementary charge [coulombs]';
C.epsilon0 = 8.85418782e-12;
C.info.epsilon0 = 'Electric constant [m-3 kg-1 s4 A2]';
end