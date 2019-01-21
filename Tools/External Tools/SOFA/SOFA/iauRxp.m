%  - - - - - - -
%   i a u R x p
%  - - - - - - -
%
%  Multiply a p-vector by an r-matrix.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     r(3,3)             r-matrix
%     p(3)               p-vector
%
%  Returned:
%     rp(3)              r * p
%
%  Note:
%     It is permissible for p and rp to be the same array.
%
%  Called:
%     iauCp        copy p-vector
%
%  This revision:  2008 October 28
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rp = iauRxp(r, p)

% Matrix r * vector p.
wrp = zeros(3,1);
for j = 1:3
    w = 0;
    for i = 1:3
        w = w + r(j,i) * p(i);
    end
    wrp(j) = w;
end

% Return the result.
rp = iauCp(wrp);

