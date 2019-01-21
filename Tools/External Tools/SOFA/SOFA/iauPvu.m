%  - - - - - - -
%   i a u P v u
%  - - - - - - -
%
%  Update a pv-vector.
%
%  This function is part of the International Astronomical Union's
%  SOFA (Standards Of Fundamental Astronomy) software collection.
%
%  Status:  vector/matrix support function.
%
%  Given:
%     dt                  time interval
%     pv(2,3)             pv-vector
%
%  Returned:
%     upv(2,3)            p updated, v unchanged
%
%  Notes:
%
%  1) "Update" means "refer the position component of the vector
%     to a new date dt time units from the existing date".
%
%  2) The time units of dt must match those of the velocity.
%
%  3) It is permissible for pv and upv to be the same array.
%
%  Called:
%     iauPpsp      p-vector plus scaled p-vector
%     iauCp        copy p-vector
%
%  This revision:  2008 November 17
%
%  SOFA release 2012-03-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function upv = iauPvu(dt, pv)

upv(1,:) = iauPpsp(pv(1,:), dt, pv(2,:));
upv(2,:) = iauCp(pv(2,:));

