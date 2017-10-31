function UT = datestr2unix(datestr,varargin)
% datastr2unix
% by John Swoboda
% This changes a datestr or a cellarray of date strings to a scalar or a
% vector of unix times.
% Inputs 
% datestr - A string or cell array of strings. 
% dateformat(Optional) - A string that shows the format of the date strings used.
% Output
% UT - A scalar or vector of unix times. 
if nargin ==1
    
    UT = 24*3600*(datenum(datestr)-datenum('jan-01-1970'));
elseif nargin==2
    UT = 24*3600*(datenum(datestr,varargin{1})-datenum('jan-01-1970'));
else
    error('Wrong number of inputs');
end