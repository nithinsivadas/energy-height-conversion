function [ error_in_q, q, qTime ] = get_error_in_q( Ne, dNe, alt, time, mode, alpha)

%calcError.m Calculates the error 

%   Input
%-----------
%   Ne     : electron density [m-3] (mean electron density)
%   dNe    : error in electron density [m-3] (standard deviation)
%   alt    : altitude [km]
%   time   : time vector [Matlab time]
%   alpha  : effective recombination rates [m^3 s^-^1]
%   mode   : Choose between the following calculations
%		'1': q = dn/dt + alpha*ne^2 (Default)
%		'2': q = alpha*ne^2 

%  Output
%-----------
%    error_in_q    : error in production rate [m-3 s-1]
%    q             : production rate [m-3 s-1]
%    qTime         : time instances where q is estimated [matlab units]

	if nargin<5
		mode = 1;
	end;

	if nargin<6
	    [q, qTime, alpha] = get_production_rate(Ne, alt, time, mode);
	else
	    [q, qTime, alpha] = get_production_rate(Ne, alt, time, mode, alpha);
	end;
	
	ialt = 1:1:length(alt);
	itime = 1:1:length(qTime);

	%% Estimating the error or standard deviation in the production rate

	% Error in the dn/dt term; But not sure about it
	ddNe2 = dNe(ialt,3:1:end).^2+dNe(ialt,1:1:end-2).^2; % Assuming 0 covariance between time instances
	dq_1  = ddNe2./(2*DT);

	% Error in alpha*Ne^2 term
	dq_2   = 3*diag(alpha.^2)*((dNe(ialt,itime).^2).*(2*ne(ialt,itime).^2+dNe(ialt,itime).^2));

	if mode == 1
		dq = (dq_1 + dq_2).^0.5
	elseif mode == 2
		dq = (dq_2).^0.5; %  Standard deviation only considering error
	end;

    error_in_q = dq;

    [isThereNAN, totalNAN] = check_nan(error_in_q);
end






