% Routine to estimate electron density along beams at
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
pfisrGDcopy = copy(pfisrGD);
pfisrGD.interpolate(pfisrGDcopy.dataloc,'Cartesian','natural');

coords = pfisrGD
%% Extracting data along beams
coordinateNo=1; % Using slant range to evaluate the number of beams
nBeams = calculate_nBeams(pfisrGD, coordinateNo);
    
for iBeam=1:1:nBeams
    [electronDensity(:,iBeam,:), xEast(:,iBeam), yNorth(:,iBeam), altitude(:,iBeam), time(:,iBeam)] =...
        extract_pfisr_beam_data(pfisrGD, 'ne', iBeam, nBeams);
    for itime=1:1:length(time(:,iBeam))
        electronDensity(:,iBeam,itime) = interp_nans(electronDensity(:,iBeam,itime));
    end;
end;

%% Inverting electron density to energy spectra
coordinate_description='[xEast, yNorth, altitude] in km';
save('2D_energy_spectra.mat','electronDensity', 'xEast', 'yNorth','altitude', 'nBeams','coordinate_description');
%% plotting
%% Need plot_2D_electron_density_slice(..) 
% figure; plot_2D_energy_slice(data, magcoords, energyBin, nBeams, 1, 110, 100);
% colormap(inferno);
% colorbar;
