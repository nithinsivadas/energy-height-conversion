%% Calculating the Minimum Energy of Scatter measured by THEMIS

%% Load Files
load('G:\My Drive\Research\Projects\Paper 2\Data\MAT Files\thmLossConeAndTrappedFlux_2008_03_26_TS96.mat');
omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,'26-Mar-2008 08:00','26-Mar-2008 13:00');

%% Calculate L-shell or equatorial parameters
smoothFactor = 10;
kext = find_irbem_magFieldModelNo('TS96');
maginput=filter_irbem_maginput(kext,maginput);
time = padData.thd.time;
xGSE = padData.thd.XYZ_GSE;
global GEOPACK1;
C = define_universal_constants;
PARMOD = zeros(10,1);
PARMODT96 = zeros(10,1);
tic
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    t = datetime(datevec(time(iTime)));
    thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
    GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));
    kp = round(thisMaginput(1)/10)+1;
    PARMOD(1) = kp;
    PARMODT01(1) = thisMaginput(5);
    PARMODT01(2) = thisMaginput(2);
    PARMODT01(3) = thisMaginput(6);
    PARMODT01(4) = thisMaginput(7);
    PARMODT01(5) = thisMaginput(8);
    PARMODT01(6) = thisMaginput(9);
    
    xGSM = onera_desp_lib_coord_trans(xGSE(iTime,:),'gse2gsm',thisTime);
    [~,~,~,X1,Y1,Z1,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
        1.0,50,(6371.2+110)/6371.2,kp,PARMOD,'T89','GEOPACK_IGRF_GSM');
    
    [~,~,~,X2,Y2,Z2,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
        -1.0,50,(6371.2+110)/6371.2,kp,PARMOD,'T89','GEOPACK_IGRF_GSM');
    
    XX{iTime} = [fliplr(X1),X2];
    YY{iTime} = [fliplr(Y1),Y2];
    ZZ{iTime} = [fliplr(Z1),Z2];
    
    Kc{iTime} = geopack_find_curvature(XX{iTime},YY{iTime},ZZ{iTime});
    
    [maxKc,indx] = max(Kc{iTime});
    minRc(iTime) = 1./maxKc;
    magEq{iTime} = [XX{iTime}(indx),YY{iTime}(indx),ZZ{iTime}(indx)];
    [BX,BY,BZ] = T89(kp,PARMOD,GEOPACK1.PSI,magEq{iTime}(1),magEq{iTime}(2),magEq{iTime}(3));
    KE(iTime) = ((((minRc(iTime).*C.RE ).*(C.e).*(BZ*10^-9)).^2).*(2^-7).*(C.me).^-1).*(10^-3).*(C.e).^-1; %keV
end
toc
%%
padData.thd.KESmooth = conv(KE,ones(smoothFactor,1)/smoothFactor,'same');
