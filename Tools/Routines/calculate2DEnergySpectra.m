% Routine to estimate energy spectra along magnetic field aligned beams at
% a particular time

%% Initializing
clear 
clc
fileName = '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
pfisrGD  = GeoData(@readMadhdf5,fileName,{'popl'});

%% Change the log electron density to linear density
antilog = @(ex,base)(base.^ex);
pfisrGD.changedata('popl','ne',antilog,{10});

%% Create magnetic field aligned data for inversion
pfisrGDmag = copy(pfisrGD);
[pfisrGDmag, magcoords, atTime, atAltitude] = geodata_magnetic_field_interpolation...
    (pfisrGDmag, '26 Mar 2008 8:00','26 Mar 2008 13:00', 60, 1);

%% Extracting data along beams
coordinateNo=1; % Using slant range to evaluate the number of beams
nBeams = get_nbeams(pfisrGD, coordinateNo);
    
for iBeam=1:1:nBeams
    [electronDensity(:,iBeam,:), xEast(:,iBeam), yNorth(:,iBeam), altitude(:,iBeam), time(:,iBeam)] =...
        get_pfisr_beam_data(pfisrGDmag, 'ne', iBeam, nBeams);
    for itime=1:1:length(time(:,iBeam))
        electronDensity(:,iBeam,itime) = interp_nans(electronDensity(:,iBeam,itime));
    end;
end;

% Note xEast, yNorth ought to be the same as magcoords(:,1) & (:,2)
%% Inverting electron density to energy spectra
energyBin = logspace(3,6,25)';
A =get_energy_dep_matrix(altitude(:,1),energyBin,...
                65,-147.5,datenum([2008 03 26 10 00 00])); %[m^-1 eV] 

            
for iBeam=1:1:nBeams
    q(:,iBeam,:) = get_production_rate(squeeze(electronDensity(:,iBeam,:)),altitude(:,1),...
        time(:,iBeam), 2);
    
    data(iBeam) = get_inverted_flux(squeeze(q(:,iBeam,:)), time(:,iBeam), altitude(:,1),...
        energyBin,A);
    
    data(iBeam).electronDensity = electronDensity(:,iBeam,:);
    
end;
coordinate_description='magcoords  - [xEast, yNorth, zUp] in km';
save('2D_energy_spectra.mat','data', 'magcoords', 'energyBin', 'nBeams','coordinate_description');
%% plotting
figure; plot_2D_energy_slice_geodetic(data, magcoords, energyBin, nBeams, 1, 110, 100);
colormap(inferno);
colorbar;
