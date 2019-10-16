%% Getting electron density and other terms for Antti Kero
% In order to estimate the chemical composition of the D-region 
clear all;

T = read_h5_data('G:\My Drive\Research\Projects\Paper 3\Data\3339_20080326_pfisrData.h5');
nLength = size(T,1);

for i = 1:1:nLength
    data1.(string(T.Name(i))) = T.Data{i};
end

data.alt = data1.alt;
data.az = data1.az;
data.el = data1.el;
data.lat = data1.lat;
data.lon = data1.lon;
data.range = data1.range;
data.time = data1.time;
data.Ne = data1.Ne;
data.q = data1.q;
data.dNeError = data1.dNeFrac;
data.energyFlux = data1.energyFlux;
data.energyBin = data1.energyBin;
data.dEnergyFluxError = data1.dEnergyFlux;
data.sensorLoc = [data1.latitude, data1.longitude, data1.altitude];

data.units.alt = '[km]';
data.units.az = '[deg]';
data.units.el = '[deg]';
data.units.lat = '[deg North]';
data.units.lon = '[deg East]';
data.units.range = '[km]';
data.units.time = '[matlab date units]';
data.units.Ne = '[m^-^3]';
data.units.q = '[m^-^3 s^-1]';
data.units.dNeError = '[m^-^3 s^-^1]';
data.units.energyFlux = '[eV/ m^2 sr s eV]';
data.units.energyBin = '[eV]';
data.units.dEnergyFluxError = '[eV/ m^2 sr s eV]';
data.units.sensorLoc = '[deg North], [deg East], [m]';

data.description.alt = 'Altitude';
data.description.az = 'Azimuth';
data.description.el = 'Elevation';
data.description.lat = 'Latitude';
data.description.lon = 'Longitude';
data.description.range = 'Slant Range along the beam';
data.description.time = 'Time in UTC';
data.description.Ne = 'Electron density';
data.description.q = 'Production rate';
data.description.dNeError = 'Error in electron density measurement';
data.description.energyFlux = 'Energy flux estimated from Max Entropy Method';
data.description.energyBin = 'Energy bin values';
data.description.dEnergyFluxError = 'Error in energy flux';
data.description.sensorLoc = 'Latitude, Longitude, and Altitude of the radar';

