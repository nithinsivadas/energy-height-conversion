function out = flatten_struct(in,n,sep)

% flatten_struct - flatten structure
% ----------------------------------
%
% out = flatten_struct(in,n,sep)
%
% Input:
% ------
%  in - arbitrary struct
%  n - number of levels to flatten (def: -1, flatten all)
%  sep - field separator (def: '__', double underscore)
%
% Output:
% -------
%  out - flat struct

%--------------------------------
% Author: Harold Figueroa
%--------------------------------
% $Revision: 1597 $
% $Date: 2005-08-17 18:33:05 -0400 (Wed, 17 Aug 2005) $
%--------------------------------

% TODO: stop recursive flattening if we hit 'namelengthmax'

% TODO: figure out how to extend to struct arrays

%-------------------------------------------------
% HANDLE INPUT
%-------------------------------------------------

%--
% set default separator
%--

if (nargin < 3)
	sep = '__';
end

%--
% set number of levels to flatten
%--

if ((nargin < 2) || isempty(n))
	n = -1;
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

%-------------------------------------------------
% FLATTEN STRUCT
%-------------------------------------------------

out = in; 

flag = 1;

switch (n)
	
	%--
	% identity
	%--
	
	case (0)	
		
	%--
	% maximally flatten
	%--
	
	case (-1)
		
		while(flag)
			[out,flag] = flatten_struct_once(out,sep);
		end
						
	%--
	% flatten fixed number of levels
	%--
	
	otherwise
			
		for k = 1:n

			[out,flag] = flatten_struct_once(out,sep); 

			if (~flag) 
				return;	
			end

		end
		
end
			

%------------------------------------------
% FLATTEN_STRUCT_ONCE
%------------------------------------------

function [out,flag] = flatten_struct_once(in,sep)

% flatten_struct_once - flatten structure one level
% -------------------------------------------------
%
% [out,flag] = flatten_struct_once(in,sep)
%
% Input:
% ------
%  in - input structure
%
% Output:
% -------
%  out - flatter structure

%--
% flatten struct one level
%--

flag = 0;

field = fieldnames(in);

for k = 1:length(field)
	
	tmp = in.(field{k});

	if (isstruct(tmp) && (length(tmp) == 1))
		
		%--
		% flag change in struct
		%--
		
		flag = 1;
		
		%--
		% get field fields and flatten field
		%--
		
		tmp_field = fieldnames(tmp);
		
		% NOTE: we string together the level fields with a double underscore

		for j = 1:length(tmp_field)
			out.([field{k}, sep, tmp_field{j}]) = tmp.(tmp_field{j});
		end
		
	%--
	% output other fields unchanged
	%--
	
	else
		
		out.(field{k}) = tmp;
		
	end
	
end
