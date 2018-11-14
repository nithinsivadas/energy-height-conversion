function [A,B,C] = GEOPACK_GSMGSE (X,Y,Z,J)
% function [A,B,C] = GEOPACK_GSMGSE (X,Y,Z,J)
%       SUBROUTINE GSMGSE (XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,J)
% C
% C CONVERTS GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDS TO SOLAR ECLIPTIC (GSE) ONES
% C   OR VICA VERSA.
% C                    J>0                J<0
% C-----INPUT: J,XGSM,YGSM,ZGSM    J,XGSE,YGSE,ZGSE
% C----OUTPUT:   XGSE,YGSE,ZGSE      XGSM,YGSM,ZGSM
% C
% C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GSMGSE IN TWO CASES:
% C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
% C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
% C
% C     LAST MODIFICATION:  MARCH 31, 2003
% C
% C     AUTHOR:  N. A. TSYGANENKO
% C

%       COMMON /GEOPACK1/ A(12),SHI,CHI,AB(13),BA(8)
% C
if J<0,
    [A,B,C] = GEOPACK_GSE2GSM(X,Y,Z);
else
    [A,B,C] = GEOPACK_GSM2GSE(X,Y,Z);
end

function [XGSE,YGSE,ZGSE] = GEOPACK_GSM2GSE(XGSM,YGSM,ZGSM)
global GEOPACK1;
XGSE=XGSM;
YGSE=YGSM*GEOPACK1.CHI-ZGSM*GEOPACK1.SHI;
ZGSE=YGSM*GEOPACK1.SHI+ZGSM*GEOPACK1.CHI;

function [XGSM,YGSM,ZGSM] = GEOPACK_GSE2GSM(XGSE,YGSE,ZGSE)
global GEOPACK1;
XGSM=XGSE;
YGSM=YGSE*GEOPACK1.CHI+ZGSE*GEOPACK1.SHI;
ZGSM=ZGSE*GEOPACK1.CHI-YGSE*GEOPACK1.SHI;

% end of function GSMGSE
% C
% C=====================================================================================
% C


