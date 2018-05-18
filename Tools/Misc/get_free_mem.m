function [mem, unit] = get_free_mem()
% usage: [mem,unit] = get_free_mem()
[~,out] = system('vmstat -s -S M | grep "free memory"');
mem = sscanf(out, '%f free memory');
unit = 'MB';
end