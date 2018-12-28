%% Estimating L-shell of 
clear all;

load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v1.h5';

omniH5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
timeMinStr = '26 Mar 2008 8:00';
timeMaxStr = '26 Mar 2008 13:00';
iBeam = 13;

[maginput,maginputTime] = generate_maginput(omniH5,timeMinStr,timeMaxStr);
maginput = interp_nans(maginput);
magFieldNo = find_irbem_magFieldModelNo('TS96');
maginput = filter_irbem_maginput(magFieldNo,maginput);

time = linspace(datenum(timeMinStr),datenum(timeMaxStr),200);
magGeodeticLatLon = h5read(fileStr,'/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt');
%% Yiqun Yu's Conductivities
baseStr = 'G:\My Drive\Research\Projects\Yiqun Yu\GLOW_output.tar\GLOW_output\GLOW_output_one_point_at_';
yiqun.Ar.tha=extract_yiqun_ascii([baseStr,'lc_tha_Ar.dat']);
yiqun.Ar.thd=extract_yiqun_ascii([baseStr,'lc_thd_Ar.dat']);
yiqun.Ar.the=extract_yiqun_ascii([baseStr,'lc_the_Ar.dat']);
yiqun.Li.tha=extract_yiqun_ascii([baseStr,'lc_tha_Li.dat']);
yiqun.Li.thd=extract_yiqun_ascii([baseStr,'lc_thd_Li.dat']);
yiqun.Li.the=extract_yiqun_ascii([baseStr,'lc_the_Li.dat']);
yiqun.Yi.tha=extract_yiqun_ascii([baseStr,'lc_tha_Yi.dat']);
yiqun.Yi.thd=extract_yiqun_ascii([baseStr,'lc_thd_Yi.dat']);
yiqun.Yi.the=extract_yiqun_ascii([baseStr,'lc_the_Yi.dat']);

%%
thisTimeStr = '26 Mar 2008 10:20';
pfisrTimeIndx = find_time(time,thisTimeStr);
% hold on; plot(datetime(datevec(yiqun.Ar.tha.time)),tha.Lm);
hold on; plot(datetime(datevec(yiqun.Ar.thd.time)),conv(thd.Lm,ones(10,1)/10,'same'),'b');
hold on; plot(datetime(datevec(yiqun.Ar.the.time)),conv(the.Lm,ones(10,1)/10,'same'),'g'); 
legend('pfisr','thd','the'); ylabel('Lm');
grid off;
title('TS96 Lm Variation with Time');
set(gca,'YScale','log');
xlim([datetime('2008-03-26 09:00') datetime('2008-03-26 12:30')]);

function [Lm,MLT,MLAT,MLON,MLT_AACGM] = get_pfisr_magnetic_coordinates(time,maginput,GDZ,magFieldNo)
    [Lm,~,~,~,~,MLT] = onera_desp_lib_make_lstar(magFieldNo,[0,0,0,0,0],0,time,GDZ(3),GDZ(1),GDZ(2),maginput);
    tup=py.aacgmv2.wrapper.get_aacgm_coord(GDZ(1), GDZ(2), GDZ(3), time, 'TRACE');
    MLAT = double(py.array.array('d',py.numpy.nditer(tup{1})));
    MLON = double(py.array.array('d',py.numpy.nditer(tup{2})));
    MLT_AACGM = double(py.array.array('d',py.numpy.nditer(tup{3})));
    Lm = abs(Lm);
end

function [Lm,MLT,MLAT,MLON,MLT_AACGM] = get_themis_magnetic_coordinates(time,maginput,GSE,magFieldNo)
    [Lm,~,~,~,~,MLT] = onera_desp_lib_make_lstar(magFieldNo,[0,0,0,0,0],3,time,GSE(1),GSE(2),GSE(3),maginput);
    tup=py.aacgmv2.wrapper.get_aacgm_coord(GSE(1), GSE(2), GSE(3), time, 'TRACE');
    MLAT = double(py.array.array('d',py.numpy.nditer(tup{1})));
    MLON = double(py.array.array('d',py.numpy.nditer(tup{2})));
    MLT_AACGM = double(py.array.array('d',py.numpy.nditer(tup{3})));
    Lm = abs(Lm);
end