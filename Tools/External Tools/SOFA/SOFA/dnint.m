% dnint(A) - round to nearest whole number (double)
function dnintA = dnint(A)

if (A<0)
    dnintA = ceil(A-0.5);
else
    dnintA = floor(A+0.5);
end