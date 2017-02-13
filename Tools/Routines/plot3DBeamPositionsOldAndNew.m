% Routine to plot existing beams positions in 3D, along with sample and
% the magnetic field aligned beam positions
%% Plotting data points along beam
clear 
clc
fileName = '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
pfisrGD  = GeoData(@readMadhdf5,fileName,{'popl'});

%% Change the log electron density to linear density
antilog = @(ex,base)(base.^ex);
pfisrGD.changedata('popl','ne',antilog,{10});

%% Converting beam data locations from AER to ENU (East North Up) coordinates
dataLocation = pfisrGD.dataloc;
az   = dataLocation(:,2);
elev = dataLocation(:,3);
slant= dataLocation(:,1);

[yNorth, xEast, zDown] = aer2ned(az,elev,slant);
 zUp = -zDown;
 
%% Creating magnetic field aligned beam locations in ENU coordinates
pfisrGDmag = copy(pfisrGD);
[pfisrGDmag, magcoords, atTime, atAltitude] = geodata_magnetic_field_interpolation...
    (pfisrGDmag, '26 Mar 2008 11:48', '26 Mar 2008 11:48', 60, 1);

x1MagEast=magcoords(:,1);
y1MagNorth=magcoords(:,2);
zMagUp = magcoords(:,3);
% Finding corresponding geodetic coordinates
pfisrloc=[65.12992 -147.47104 213];
[maglat, maglon, magh] = ned2geodetic...
    (magcoords(:,2), magcoords(:,1), -magcoords(:,3),pfisrloc(1),pfisrloc(2),pfisrloc(3),wgs84Ellipsoid('km'));
 
 %% Time reduce
pfisrGD.timereduce(344);

%% Interpolation of data along a 3D grid in order to feasibly plot
xvec = linspace(-100,200,50);
yvec = linspace(-100,200,50);
zvec = linspace(50,200,50);

[X,Y,Z] = meshgrid(xvec,yvec,zvec);

[Xmat,Ymat] = meshgrid(xvec,yvec);
newcoords = [X(:),Y(:),Z(:)];
pfisrGDorig = copy(pfisrGD);

pfisrGD.interpolate(newcoords,'Cartesian','natural');

%% Plot Data

k=1;
irt =1;
curirt=1;
xlist = [0];
ylist = [0];
zlist = [60,110];

% make a four quadrant figure with the combined data in 3-D, two
% 2-D altitude slices and a beam pattern.
hfig = figure('Color',[1,1,1],'Position',[680,250,1100,725]);

%% Figure 1
hax = subplot(1,1,1);
pfisrslic = sliceGD(pfisrGD,xlist,ylist,zlist,'key','ne','Fig',hfig,'axh',hax,'title'...
    ,'','time',curirt,'bounds',[5e9,5e11]);
hcb = colorbar();
ylabel(hcb,'N_e in m^{-3}');
 axis tight
 view(-40,30);

%%
hold on;
scatter3(xEast, yNorth, zUp,'.');
xlabel('East [km]');
ylabel('North [km]');
zlabel('Altitude [km]');

hold on;
scatter3(x1MagEast, y1MagNorth, zMagUp,'.r');

%% Using Magnetic Field aligned pfisrGD object to extract data along a beam
[msrParValue, altitudeGrid, time] = get_pfisr_beam_data(pfisrGDmag, 'ne', 1);
