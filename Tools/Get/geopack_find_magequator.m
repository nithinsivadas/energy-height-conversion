function [magEqGEO] = geopack_find_magequator(magFieldNo,...
    maxLength,sysaxes,thisTime,x1,x2,x3,maginput,runGEOPACK_RECALC)
%geopack_find_magequator Find magnetic equator points using GEOPACK tracing

if nargin<9
    runGEOPACK_RECALC = true;
end

if sysaxes ~=2
    xGSM = onera_desp_lib_coord_axes([x1,x2,x3],[sysaxes 2],thisTime);
else
    xGSM = [x1,x2,x3];
end

if isempty(maxLength)
    maxLength = 50;
end

[PARMOD,IOPT,EXNAME] = get_parmod(magFieldNo,maginput);
DIR = get_dir(xGSM,hemiflag);
RLIM = maxLength;

INNAME = 'GEOPACK_IGRF_GSM';

if runGEOPACK_RECALC == true
    t = datetime(thisTime,'ConvertFrom','datenum');
    GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));
end

[~,~,~,XX,YY,ZZ,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
    DIR,RLIM,1,IOPT,PARMOD,EXNAME,INNAME);

[Kc] = geopack_find_curvature(XX,YY,ZZ);

[~,indx] = max(Kc);

magEqGSM = [XX(indx),YY(indx),ZZ(indx)];
magEqGEO = onera_desp_lib_coord_axes(magEqGSM,[2 sysaxes],thisTime);

end

function [PARMOD,IOPT,magStr] = get_parmod(magFieldNo,maginput)

PARMOD = zeros(10,1);

if magFieldNo==4
    kp = round(maginput(1)/10)+1;
    magStr = 'T89';
    PARMOD(1) = kp;
    IOPT = kp;
elseif magFieldNo==7
    magStr = 'T96';
    PARMOD(1) = maginput(5);
    PARMOD(2) = maginput(2);
    PARMOD(3) = maginput(6);
    PARMOD(4) = maginput(7);
    IOPT = 0;
elseif magFieldNo==9
    magStr = 'T01';
    PARMOD(1) = maginput(5);
    PARMOD(2) = maginput(2);
    PARMOD(3) = maginput(6);
    PARMOD(4) = maginput(7);
    PARMOD(5) = maginput(8);
    PARMOD(6) = maginput(9);
    IOPT=0;
else
    error('GEOPACK version does not contain the specified external field.');
end

end

function DIR = get_dir(xGSM,hemiflag)

DIR = sign(xGSM(3));
    if hemiflag == 0
        DIR = DIR*-1;
    elseif hemiflag == +1
        DIR = -1*abs(DIR);
    elseif hemiflag == -1
        DIR = abs(DIR);
    elseif hemiflag == 2
    else
        error('hemiflag value out of range');
    end

end
