function B = interp_nans(A)
% This function removes nan by interpolating along altitude
%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------

	x=1:1:size(A,1);

	for ty=1:1:size(A,2)
	
	    y=A(:,ty);
	    xi=x(find(~isnan(y)));
	    yi=y(find(~isnan(y)));
	    B(:,ty)=interp1(xi,yi,x,'linear','extrap');
	
	end;

	[isThereNAN, totalNAN] = check_nan(B);

end