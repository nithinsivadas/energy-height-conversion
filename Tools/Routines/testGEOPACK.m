%% Test GEOPACK & IRBEM
% In conclusion GEOPACK works reasonable for TS89
clear all;
%% Initializing
poesMatFile = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326.mat';
load(poesMatFile);
omniHDF5File = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
[maginput,timeMaginput]=generate_maginput(omniHDF5File,'26 Mar 2008 11:00','26 Mar 2008 12:00');
%%
timeStr = '26-Mar-2008 11:00:00';
t = datetime(timeStr,'InputFormat','dd-MMM-yyyy HH:mm:ss');

GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));

thisTime = datenum(timeStr);
alt = 850;
lat = interp1(poes.time,poes.lat,thisTime);
lon = interp1(poes.time,poes.lon,thisTime);
xGEO = onera_desp_lib_rotate([alt, lat, lon],'gdz2geo');
xGSM = onera_desp_lib_rotate(xGEO,'geo2gsm',thisTime);
thisMaginput = interp1(timeMaginput',maginput,thisTime);
kp = round(thisMaginput(1)/10)+1;
PARMOD = zeros(10,1);
PARMOD(1) = kp;
PARMODT96 = zeros(10,1);
PARMODT96(1) = thisMaginput(5);
PARMODT96(2) = thisMaginput(2);
PARMODT96(3) = thisMaginput(6);
PARMODT96(4) = thisMaginput(7);

[XF,YF,ZF,XX,YY,ZZ,L] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
    -1.0,50,(6371.2+103)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');

xFGEO = onera_desp_lib_rotate([XF,YF,ZF],'gsm2geo',thisTime);
xFGDZ = onera_desp_lib_rotate(xFGEO,'geo2gdz');

NFoot = onera_desp_lib_find_foot_point...
            (find_irbem_magFieldModelNo('TS96'),[0,0,0,0,0],2,thisTime,...
            xGSM(1),...
            xGSM(2),...
            xGSM(3),...
            110,+1,thisMaginput);
        
fprintf('GEOPACK: %6.2f deg %6.2f deg %6.2f km \n',xFGDZ(2),xFGDZ(3),xFGDZ(1));
fprintf('IRBEM  : %6.2f deg %6.2f deg %6.2f km \n',NFoot(2),NFoot(3),NFoot(1));
