function output = struct_trim(input,data_keep)
% struct_trim.m
% by John Swoboda
% This function will take a flat struct as input and go throught he arrays
% and trim them down according to the logical array or arrya of indices
% data_keep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%
% input - A struct made up of arrays of the same size.
% data_keep - A logical array of the same size of the arrays in input or an
% array of indices with the max being the number of elements in the arrays
% in input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% 
% output - A struct made up of arrays of the same sizes with the same field
% names as in input, but with less data.

fields = fieldnames(input);
output = struct();
for k = 1:length(fields);
    tmp = input.(fields{k});
    output.(fields{k}) = tmp(data_keep);
end
