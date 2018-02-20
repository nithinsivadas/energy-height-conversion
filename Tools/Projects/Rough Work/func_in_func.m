function [result] = func_in_func(varargin)
%FUNC_IN_FUNC Summary of this function goes here
%   see: http://stackoverflow.com/questions/2730029/
x = varargin{1};
y = varargin{2};
feq = varargin{3};
result = feq(x,y);
end

