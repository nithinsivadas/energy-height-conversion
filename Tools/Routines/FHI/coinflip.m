
function y=coinflip(t,p,repeat)

% Running an experiment where events occur according to an underlying uniformly 
% distributied random number, such that the probability of occurence is $p$. 
%
% Input: 
% t      - is the time
% p      - is the probability that an event occurs at this time instant
% repeat - repeat the experiment 'repeat' number of times
%          [default: 1]
% Output: 
% y       - is a [TxRepeat] matrix showing the status of the event (success/failure) 
%           at each time instant (row), for each experiment (column)

    if nargin < 3
        repeat = 1;
    end
        
    n = length(t);
    y = zeros(n,repeat);
    x = rand(n,repeat);
    y(x<p) = 1;
end