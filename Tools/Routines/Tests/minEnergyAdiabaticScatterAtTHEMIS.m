%% Calculating the Minimum Energy of Scatter measured by THEMIS

%% Load Files
tic
disp('Loading...');
load('G:\My Drive\Research\Projects\Paper 2\Data\MAT Files\thmLossConeAndTrappedFlux_2008_03_26_TS96.mat');
omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,'26-Mar-2008 08:00','26-Mar-2008 13:30');
toc
%% Calculate L-shell or equatorial parameters
tic
disp('Filtering...');
smoothFactor = 10;
kext = find_irbem_magFieldModelNo('TS96');
maginput = interp_nans(maginput);
maginput=filter_irbem_maginput(kext,maginput);
toc
%% THD
tic
disp('THD...');
time = padData.thd.time;
xGSE = padData.thd.XYZ_GSE;
[KE_thd,minRc_thd]=get_min_scatter_energy(time,magTime,maginput,xGSE);
toc
%% THE
tic
disp('THE...');
time = padData.the.time;
xGSE = padData.the.XYZ_GSE;
[KE_the,minRc_the]=get_min_scatter_energy(time,magTime,maginput,xGSE);
toc
%%
padData.thd.KE = KE_thd;
padData.thd.KESmooth = conv(KE_thd,ones(smoothFactor,1)/smoothFactor,'same');

padData.the.KE = KE_the;
padData.the.KESmooth = conv(KE_the,ones(smoothFactor,1)/smoothFactor,'same');

%% Plot
disp('Plotting...');
timeMinStr = '26-Mar-2008 08:30';
timeMaxStr = '26-Mar-2008 12:30';
totalPanelNo=1;
p = create_panels(figure,'totalPanelNo',totalPanelNo,'margintop',4,'panelHeight',30);
q=p(1);
q(1).select();
scatter(padData.thd.time,padData.thd.KE,10,'filled','MarkerFaceColor','b', 'MarkerFaceAlpha',.4);
% fthd=fit(padData.thd.time,padData.thd.KE','smoothingspline','SmoothingParam',0.99991);
% hold on;
% ax1=plot(padData.thd.time,fthd(padData.thd.time),'b','LineWidth',1);
hold on;
ax1=plot(padData.thd.time,movmean(padData.thd.KE,10),'b','LineWidth',1);

% fthe=fit(padData.the.time,padData.the.KE','smoothingspline','SmoothingParam',0.999999);
scatter(padData.the.time,padData.the.KE,10,'filled','MarkerFaceColor','c', 'MarkerFaceAlpha',.4);
% ax2=plot(padData.the.time,fthe(padData.the.time),'c','LineWidth',1);
hold on;
ax2=plot(padData.the.time,movmean(padData.the.KE,10),'c','LineWidth',1);

set(gca,'YScale','log','YTick',[1,10,100,300],...
    'YGrid','on','YMinorGrid','off','YLim',[1 1000]);
ylabel({'min e^- energy','scattered','[keV]'});
legend([ax1, ax2],{'THEMIS-D','THEMIS-E'},'Location','southwest');
label_time_axis(padData.the.time,true,0.5,timeMinStr,timeMaxStr);

%% Function 
function [KE,minRc]=get_min_scatter_energy(time,magTime,maginput,xGSE)
    global GEOPACK1;
    C = define_universal_constants;
    PARMOD = zeros(10,1);
    PARMODT96 = zeros(10,1);
    for iTime = 1:1:length(time)
        thisTime = time(iTime);
        t = datetime(datevec(time(iTime)));
        thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
        GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));
        PARMODT96(1) = thisMaginput(5);
        PARMODT96(2) = thisMaginput(2);
        PARMODT96(3) = thisMaginput(6);
        PARMODT96(4) = thisMaginput(7);

        xGSM = onera_desp_lib_coord_trans(xGSE(iTime,:),'gse2gsm',thisTime);
        [~,~,~,X1,Y1,Z1,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
            1.0,50,(6371.2+110)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');

        [~,~,~,X2,Y2,Z2,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
            -1.0,50,(6371.2+110)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');

        XX{iTime} = [fliplr(X1),X2];
        YY{iTime} = [fliplr(Y1),Y2];
        ZZ{iTime} = [fliplr(Z1),Z2];

        Kc{iTime} = geopack_find_curvature(XX{iTime},YY{iTime},ZZ{iTime});

        [maxKc,indx] = max(Kc{iTime});
        minRc(iTime) = 1./maxKc;
        magEq{iTime} = [XX{iTime}(indx),YY{iTime}(indx),ZZ{iTime}(indx)];
        [BX,BY,BZ] = T96(0,PARMODT96,GEOPACK1.PSI,magEq{iTime}(1),magEq{iTime}(2),magEq{iTime}(3));
        KE(iTime) = ((((minRc(iTime).*C.RE ).*(C.e).*(BZ*10^-9)).^2).*(2^-7).*(C.me).^-1).*(10^-3).*(C.e).^-1; %keV
    end
end
