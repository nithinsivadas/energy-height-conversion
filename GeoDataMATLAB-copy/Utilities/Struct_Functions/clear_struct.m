function out = clear_struct(in)

% clear_struct - clear contents of struct
% ---------------------------------------
%
% out = clear_struct(in)
%
% Input:
% ------
%  in - structure
%
% Output:
% -------
%  out - structure with empty values

%--------------------------------
% Author: Harold Figueroa
%--------------------------------
% $Revision: 1600 $
% $Date: 2005-08-18 17:41:06 -0400 (Thu, 18 Aug 2005) $
%--------------------------------

%-----------------------------------
% HANDLE INPUT
%-----------------------------------

%--
% check for scalar struct input
%--

if (~isstruct(in) || (length(in) > 1))
	error('Scalar structure input is required.');
end

% NOTE: return same on empty, although this is a different type of empty

if (length(in) == 0)
	out = in; return;
end

%-----------------------------------
% CLEAR STRUCT
%-----------------------------------

%--
% get fields from flattened struct and put together input for struct
%--

in = fieldnames(flatten_struct(in));

in(:,2) = {[]};

% NOTE: we transpose to interleave fields and values in comma-separated list

in = in';

%--
% build empty struct and unflatten
%--

out = unflatten_struct(struct(in{:}));