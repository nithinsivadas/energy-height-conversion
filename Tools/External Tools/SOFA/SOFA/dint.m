% dint(A) - truncate to nearest whole number towards zero (double)
function dintA = dint(A) 

if (A<0)
    dintA = ceil(A);
else
    dintA = floor(A);
end



