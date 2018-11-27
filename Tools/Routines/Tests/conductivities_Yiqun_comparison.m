%% Routine to compare Yiqun's Run with PFISR
clear all
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v1.h5';
data = read_h5_data(fileStr,'/conductivity');
load('G:\My Drive\Research\Projects\Yiqun Yu\Data\Final\TS96\lossConeFluxTHMADE_20080326_TS96_with_footpoints.mat');
%% Parsing variables
Ne = data.Data{2};
Ne(Ne<10^6) = 10^6;

sigmaH = data.Data{3};
sigmaP = data.Data{4};
alt = data.Data{1};
time = data.Data{5};

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

%% Generating Lm, MLT,MLAT, MLON for PFISR
iBeam = 13;
timeMinStr = '26-Mar-2008 08:00';
timeMaxStr = '26-Mar-2008 13:00';
pfisrH5 = fileStr;
omniH5 = 'G:\My Drive\Research\Projects\Data\omni.h5';
magGeodeticLatLon = h5read(pfisrH5,'/magneticFieldAlignedCoordinates/magGeodeticLatLonAlt');
GDZ = interp1(magGeodeticLatLon(:,iBeam,3),squeeze(magGeodeticLatLon(:,iBeam,:)),85); %85 km projection
magFieldStr = 'TS89';
magFieldNo = find_irbem_magFieldModelNo(magFieldStr);
[maginput,maginputTime] = generate_maginput(omniH5,timeMinStr,timeMaxStr); 
maginput = interp_nans(maginput);
maginput = filter_irbem_maginput(magFieldNo, maginput);
%% PFISR
for i=1:1:length(time)
    thisTime = time(i);
    thisMaginput = interp1(maginputTime,maginput,thisTime);
    [Lm(i), MLT(i), MLAT(i), MLON(i), MLT_AACGM(i)] = ...
        get_pfisr_magnetic_coordinates(thisTime,thisMaginput,GDZ,magFieldNo);
end
%% THD
for i=1:1:length(yiqun.Ar.thd.time)
    thisTime = yiqun.Ar.thd.time(i);
    GSE = padData.thd.XYZ_GSE(i,:);
    thisMaginput = interp1(maginputTime,maginput,thisTime);
    [thd.Lm(i), thd.MLT(i), thd.MLAT(i), thd.MLON(i), thd.MLT_AACGM(i)] = ...
        get_themis_magnetic_coordinates(thisTime,thisMaginput,GSE,magFieldNo);
    
end
% THA
for i=1:1:length(yiqun.Ar.tha.time)
    thisTime =  yiqun.Ar.tha.time(i);
    GSE = padData.thd.XYZ_GSE(i,:);
    thisMaginput = interp1(maginputTime,maginput,thisTime);
    [tha.Lm(i), tha.MLT(i), tha.MLAT(i), tha.MLON(i), tha.MLT_AACGM(i)] = ...
        get_themis_magnetic_coordinates(thisTime,thisMaginput,GSE,magFieldNo);
    
end

% THE
for i=1:1:length(yiqun.Ar.the.time)
    thisTime =  yiqun.Ar.the.time(i);
    GSE = padData.the.XYZ_GSE(i,:);
    thisMaginput = interp1(maginputTime,maginput,thisTime);
    [the.Lm(i), the.MLT(i), the.MLAT(i), the.MLON(i), the.MLT_AACGM(i)] = ...
        get_themis_magnetic_coordinates(thisTime,thisMaginput,GSE,magFieldNo);
    
end
%% Plotting Ne
iBeam = 13;
% DregionAltIndx = find_altitude(alt(iBeam,:)',70):1:find_altitude(alt(iBeam,:)',107); 
% EregionAltIndx = find_altitude(alt(iBeam,:)',107)+1:1:find_altitude(alt(iBeam,:)',160);

h = figure;
timeMinStr = '26-Mar-2008 08:00';
timeMaxStr = '26-Mar-2008 12:00';
p = create_panels(h,'totalPanelNo',4,'panelHeight',30,'demargin',20);
dtime = 0.5;
% PFISR
colormap(inferno);
yaxisLim = [70 160];
cLim = [9 12];
p(1,1).select();
inputNe = squeeze(Ne(:,iBeam,:));
inputNe = interp_nans(inputNe);
plot_2D_time_series(time,alt(iBeam,:),log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb1=colorbar('eastoutside');
cb1.Position(1) = 0.9;
cb1.Position(3) = 0.01;
cb1.Label.String = '[m^-^3]';
ylim(yaxisLim);
ylabel({'PFISR N_e','Alt [km]'});
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);
[TTick,TTickLim]=label_time_axis(time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,MLAT,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,time,MLT,'MLT',2);
add_horizontal_axes(TTick,TTickLim,time,Lm,[magFieldStr,' Lm'],3);

p(1,3).select();
inputNe = yiqun.Ar.thd.NeOut * 10^6;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.thd.time,yiqun.Ar.thd.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[m^-^3]';
ylim(yaxisLim);
ylabel({'ThmD N_e','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.thd.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,thd.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,2).select();
inputNe = yiqun.Ar.tha.NeOut * 10^6;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.tha.time,yiqun.Ar.tha.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[m^-^3]';
ylim(yaxisLim);
ylabel({'ThmA N_e','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.tha.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,tha.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,4).select();
inputNe = yiqun.Ar.the.NeOut * 10^6;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.the.time,yiqun.Ar.the.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[m^-^3]';
ylim(yaxisLim);
ylabel({'ThmE N_e','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.the.time,1,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlat,'MLAT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlt,'MLT',3);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,the.Lm,[magFieldStr,' Lm'],4);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

%% Plot Conductivity SigmaH
% PFISR

h = figure;
p = create_panels(h,'totalPanelNo',4,'panelHeight',30,'demargin',20);

colormap(inferno);
yaxisLim = [70 160];
cLim = [-6 -2];
p(1,1).select();
inputNe = real(squeeze(sigmaH(:,iBeam,:)));
inputNe(inputNe<0) = 10^-8;
inputNe = interp_nans(inputNe);
plot_2D_time_series(time,alt(iBeam,:),log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb1=colorbar('eastoutside');
cb1.Position(1) = 0.9;
cb1.Position(3) = 0.01;
cb1.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'PFISR \sigma_H','Alt [km]'});
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);
label_time_axis(time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,MLAT,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,time,MLT,'MLT',2);
add_horizontal_axes(TTick,TTickLim,time,Lm,[magFieldStr,' Lm'],3);

p(1,3).select();
inputNe = yiqun.Ar.thd.sigmaH ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.thd.time,yiqun.Ar.thd.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmD \sigma_H','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.thd.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,thd.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,2).select();
inputNe = yiqun.Ar.tha.sigmaH ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.tha.time,yiqun.Ar.tha.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmA \sigma_H','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.tha.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,tha.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,4).select();
inputNe = yiqun.Ar.the.sigmaH ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.the.time,yiqun.Ar.the.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmE \sigma_H','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.the.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,the.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

%% Plot Conductivity SigmaP
% PFISR

h = figure;
p = create_panels(h,'totalPanelNo',4,'panelHeight',30,'demargin',20);

colormap(inferno);
yaxisLim = [70 160];
cLim = [-7 -3];
p(1,1).select();
inputNe = real(squeeze(sigmaP(:,iBeam,:)));
inputNe(inputNe<0) = 10^-8;
inputNe = interp_nans(inputNe);
plot_2D_time_series(time,alt(iBeam,:),log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
caxis(cLim);
cb1=colorbar('eastoutside');
cb1.Position(1) = 0.9;
cb1.Position(3) = 0.01;
cb1.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'PFISR \sigma_P','Alt [km]'});
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);
[TTick,TTickLim]=label_time_axis(time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,time,MLAT,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,time,MLT,'MLT',2);
add_horizontal_axes(TTick,TTickLim,time,Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,3).select();
inputNe = yiqun.Ar.thd.sigmaP ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.thd.time,yiqun.Ar.thd.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
colorbar;
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmD \sigma_P','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.thd.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,yiqun.Ar.thd.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.thd.time,thd.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,2).select();
inputNe = yiqun.Ar.tha.sigmaP ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.tha.time,yiqun.Ar.tha.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
colorbar;
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmA \sigma_P','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.tha.time,0,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlat,'MLAT',1);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,yiqun.Ar.tha.mlt,'MLT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.tha.time,tha.Lm,[magFieldStr,' Lm'],3);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);

p(1,4).select();
inputNe = yiqun.Ar.the.sigmaP ;
inputNe = interp_nans(inputNe);
plot_2D_time_series(yiqun.Ar.the.time,yiqun.Ar.the.alt,log10(inputNe'),dtime,0,timeMinStr,timeMaxStr); 
colorbar;
caxis(cLim);
cb2=colorbar('eastoutside');
cb2.Position(1) = 0.9;
cb2.Position(3) = 0.01;
cb2.Label.String = '[S/m]';
ylim(yaxisLim);
ylabel({'ThmE \sigma_P','Alt [km]'});
[TTick,TTickLim]=label_time_axis(yiqun.Ar.the.time,1,dtime,timeMinStr,timeMaxStr);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlat,'MLAT',2);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,yiqun.Ar.the.mlt,'MLT',3);
add_horizontal_axes(TTick,TTickLim,yiqun.Ar.the.time,the.Lm,[magFieldStr,' Lm'],4);
xlim([datenum(timeMinStr) datenum(timeMaxStr)]);


%% Trial


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