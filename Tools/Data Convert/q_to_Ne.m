function [Ne] = q_to_Ne(q, alt, time, alpha)

%% q_to_Ne.m Estimates electron density from given production rates 
%            using static recombination coefficients

% Input:
% q     : Production rate [m^-3 s^-1]
% alt   : height [length(height) x 1] [km]
% time  : time  [1 x length(time)] [s]
% alpha : User defined effective recombination coefficients 
%         per altitude point [m^3 s^-1]
%         Default: Vickerey et al., 1982

% Output:
% Ne    : Electron density [m^-3]

%%
%----------------------------------------------------------------------------
% Modified: 22nd Sep 2016 
% Created : 22nd Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%

	if nargin<4
	   alpha   = get_eff_recomb_rate(alt);
    end;
	
	Ne = (diag(alpha.^-1)*q).^0.5;

	[isThereNAN, totalNAN] = check_nan(Ne);

end