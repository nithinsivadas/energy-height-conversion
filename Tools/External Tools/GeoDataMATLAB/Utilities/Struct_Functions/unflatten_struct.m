function out = unflatten_struct(in,sep)

% unflatten_struct - unflatten struture
% -------------------------------------
%
% out = unflatten_struct(in)
%
% Input:
% ------
%  in - flattened struct
%  sep - field separator (def: '__', double underscore)
%
% Output:
% -------
%  out - unflattened struct

%--------------------------------
% Author: Harold Figueroa
%--------------------------------
% $Revision: 1597 $
% $Date: 2005-08-17 18:33:05 -0400 (Wed, 17 Aug 2005) $
%--------------------------------

%---------------------------------------
% HANDLE INPUT
%---------------------------------------

%--
% set separator
%--

if (nargin < 2) 
	sep = '__';
end

%--
% check for struct input
%--

if (~isstruct(in))
	error('Struct input is required.');
end

%--
% check for scalar struct
%--

if (length(in) ~= 1)
	error('Scalar struct input is required.');
end

%---------------------------------------
% UNFLATTEN STRUCTURE
%---------------------------------------

% NOTE: the order of flattening and unflattening matters

field = fieldnames(in);

for k = 1:length(field)
	eval(['out.', strrep(field{k},sep,'.'), ' = in.', field{k}, ';']);
end