function [A,B,C] = GEOPACK_GEIGEO(X,Y,Z,J)
% function [A,B,C] = GEOPACK_GEIGEO(X,Y,Z,J)
%       SUBROUTINE GEIGEO (XGEI,YGEI,ZGEI,XGEO,YGEO,ZGEO,J)
% C
% C   CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
% C   OR VICA VERSA.
% C                    J>0                J<0
% C----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
% C----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI
% C
% C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEIGEO IN TWO CASES:
% C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
% C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
% C
% C     LAST MODIFICATION:  MARCH 31, 2003
% 
% C     AUTHOR:  N. A. TSYGANENKO
% C

%      COMMON /GEOPACK1/ A(27),CGST,SGST,B(6)
if J<0,
    [A,B,C] = GEOPACK_GEO2GEI(X,Y,Z);
else
    [A,B,C] = GEOPACK_GEI2GEO(X,Y,Z);
end

function [XGEO,YGEO,ZGEO] = GEOPACK_GEI2GEO(XGEI,YGEI,ZGEI);
global GEOPACK1
% C
XGEO=XGEI*GEOPACK1.CGST+YGEI*GEOPACK1.SGST;
YGEO=YGEI*GEOPACK1.CGST-XGEI*GEOPACK1.SGST;
ZGEO=ZGEI;

function [XGEI,YGEI,ZGEI] = GEOPACK_GEO2GEI(XGEO,YGEO,ZGEO)
global GEOPACK1
XGEI=XGEO*GEOPACK1.CGST-YGEO*GEOPACK1.SGST;
YGEI=YGEO*GEOPACK1.CGST+XGEO*GEOPACK1.SGST;
ZGEI=ZGEO;
% end of function GEIGEO
% C
% C=======================================================================================
% C


