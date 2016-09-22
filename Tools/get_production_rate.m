function [ q, qTime, alpha ] = get_production_rate(Ne, alt, time, mode, alpha)

%% get_production_rate: Estimates the production rate from measured electron density
%  assuming a effective recombination coefficient that is static with time

% Input:
% Ne : Electron density [length(height) x length(time)] [m-3]
% alt             : height vector [length(height) x 1] [km]
% time            : time vector [1 x length(time)] [s]
% mode            : Choose between the following calculations
%		       '1': q = dn/dt + alpha*ne^2 (Default)
%		       '2': q = alpha*ne^2 
% alpha           : User defined effective recombination coefficients 
%                   per altitude point [m^3 s^-1]

% Output:
% q     		  : production rate [m-3 s-1]
% qTime           : time instances where q is estimated [s]
% alpha           : Effective recombination coefficients 
%                   used for the calculation [m^3 s^-1]

%%
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     : Semeter & Kamalabadi 2005
%----------------------------------------------------------------------------
%%

	if nargin<5
	   alpha   = eff_recombination_rate(alt);
    end;
    
    if nargin<4
		mode = 1;
    elseif nargin<3
		error('TooFewInput: At least 3 input arguments required')
	end;


	ialt = 1:1:length(alt);

	switch mode
		case 1
			timePrev = time(1:1:end-2);
			timePost = time(3:1:end);

			dt = 86400*(timePrev-timePost);
			dn = Ne(ialt,3:1:end)-Ne(ialt,1:1:end-2);
			[DT,DNE]= meshgrid(dt,Ne(ialt,1));

			dn_dt = dn./(2*DT);
			itime = 1:1:size(timePrev);
			q_1 = dn_dt(:,itime);
			q_2 = diag(alpha)*(Ne(ialt,itime).^2);
			q = q_1 + q_2; 
			qTime=timePrev;

		case 2
			itime = 1:1:size(time);
			q = diag(alpha)*(Ne(ialt,itime).^2);
			qTime = time;

		otherwise
			error('Mode not specified');
	end


end
