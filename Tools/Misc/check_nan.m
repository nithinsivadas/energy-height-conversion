function [isThereNAN, totalNAN] = check_nan(thisArray)
%% check_nan.m Check if there is any nan values in the given array
%--------------------------------------------------------------------------
% Input
%------
% thisArray : An array or a matrix
%--------------------------------------------------------------------------
% Output
%-------
% isThereNAN : A boolean variable, True => there is nan values in the
%              array, and False => there is no nan values in the array
% totalNAN   : The total number of NAN values
%
%--------------------------------------------------------------------------
% Modified: 22nd Sep 2016 
% Created : 22nd Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%--------------------------------------------------------------------------

thisArrayNAN = isnan(thisArray);
totalNAN = sum(thisArrayNAN(:));

    if totalNAN>0
        isThereNAN = true;
    else
        isThereNAN = false;
    end

    if isThereNAN==true
        warning(['There are ',num2str(totalNAN),...
        ' nan values in your data. The results may not be accurate.']);
    end
    
end

