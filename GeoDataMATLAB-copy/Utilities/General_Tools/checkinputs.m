function checkinputs(labels,in_vals,poss_labels,vals,varnames)
% checkinputs.m
% by John Swoboda 3/20/2014
% This function can be used to take a set of inputs and labels and
% determine what the variables names will be.  The variables will then be
% sent to the workspace of the calling function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% labels - A P length cell of strings, where P is the number of input
%           parameters that were used for the particular function call.
% in_vals - A P length cell that holds the values that are input.
% poss_labels - A N length cell, N is the number of possible inputs, that
%           holds the label names.
% vals - A N length cell that holds default values for the variables.
% varnames - A N length cell of strings that hold the variable names to be
%           used in the new workspace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine which variables are used.
for k =1:length(labels)
    logout = strcmpi(labels{k},poss_labels);
    if ~any(logout)
        error(['Input string ',labels{k}, ' is not a valid input label']);
    elseif sum(logout)>1
        error(['Input string ',labels{k}, ' is present more then once']);
    end
    vals{logout} = in_vals{k};
end
% send to the workspace of the calling function.
for k = 1:length(poss_labels)
    assignin('caller',varnames{k},vals{k});
end
