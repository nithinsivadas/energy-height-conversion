% dsign(A,B) - magnitude of A with sign of B (double)
function dsignAB = dsign(A,B)

if(B<0)
    dsignAB = -abs(A);
else
    dsignAB = abs(A);
end

