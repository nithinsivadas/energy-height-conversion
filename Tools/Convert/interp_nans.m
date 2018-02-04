function B = interp_nans(A)
%% interp_nans.m This function removes nan by interpolating along altitude
%--------------------------------------------------------------------------
% Input
%------
% A - Input altitude vs. time matrix [nh x nT]
%--------------------------------------------------------------------------
% Output
%-------
% B - Interpolated altitude vs. time matrix, along the altitude directon
%     with nans removed [nh x nT]
%--------------------------------------------------------------------------
% Modified: 17th Jan 2018 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
     
	x=1:1:size(A,1);
    
	for ty=1:1:size(A,2)
	
	    y=A(:,ty);
	    xi=x(find(~isnan(y)));
	    yi=y(find(~isnan(y)));
	   
        if sum(isnan(y))< size(A,1)-2
        B(:,ty)=interp1(xi,yi,x,'linear','extrap');
        else
        B(:,ty)=A(:,ty);
        end
	
    end

	[isThereNAN, totalNAN] = check_nan(B);

end