function [alpha] = eff_recombination_rate(altitudeGrid)
%% eff_recombination_coeff.m Calculates the effective recombination coefficient
% for the given altitude grid based on Vickrey et al., 1982 

% Input
% altitudeGrid : An altitude vector along which you need the effective recombination 
%                rates to be calculation 
%                [km]

% Output
% alpha        : Effective recombination rate [m^3 s^-1]

%%
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     : Vickrey et al., 1982
%----------------------------------------------------------------------------
%%

	alpha = 2.5*10^-12*exp(-altitudeGrid./51.2); % [m^3/s] %z is in [km]

end

