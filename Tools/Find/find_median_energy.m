function [ medianIndex ] = find_median_energy( array, dim )
%% find_median_energy.m Calculates the energy bin that cuts the cumulative energy flux distribution into half
%-------------------------------------------------------------------------
% Input:
%--------
% array : A 2-D array - [Energy x Time] of normalized cimulative energy
%         flux
% dim   : The energy dimension - by default 1 as per above definition of
%         array. It can also be 2
%-------------------------------------------------------------------------
% Output :
%---------
% medianIndex : The index number of the energy bin array which cuts the 
%               cumulative energy flux distribution in half
%%
%------------------------------------------------------------------------
% Modified: 7th Feb 2017 
% Created : 21st Dec 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%
  if nargin==1
    dim = 1; 
  end
  
  if dim == 2
      array=array';
  end;
  
  if dim>2
      error('Dimesion cannot be greater than 2');
  end;
  
  for thisCol=1:1:size(array,2)
      array(:,thisCol)=100*array(:,thisCol)/max(array(:,thisCol));
      medianIndex(thisCol) = find(array(:,thisCol)>50,1);
  end;

end

