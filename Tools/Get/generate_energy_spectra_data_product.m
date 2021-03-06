function status = generate_energy_spectra_data_product...
    (amisrFileStr,amisrRootPath,dascRootPath,...
    projectionAltPFISR,projectionAltDASC,energyBin,altLim,outputH5FileStr,...
    timeMinStr,timeMaxStr,...
    dascMinElevation,dascCalFileLat,dascCalFileLon, dascSetDownloadFlag)
%UNTITLED2 Generates the energy spectra data product that contains the
%magnetically field aligned electron densities, energy flux of primary
%electrons that generate the aurora, and DASC camera data to put them in
%context. 
%   Detailed explanation goes here
if nargin<13
    setDownloadDASCFlag = true;
end
if nargin<12 || isempty(dascCalFileLon)
    dascCalFileLon = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glon.dat'];
if nargin<11 || isempty(dascCalFileLat)
    dascCalFileLat = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'Toshi_Cal',filesep,'200803261040uafdasc_pkr_glat.dat'];
end
if nargin < 10
    dascMinElevation = 30;
end
if nargin < 9 || isempty(timeMaxStr)
    timeMaxStr = [];
end
if nargin < 8 || isempty(timeMinStr)
    timeMinStr = [];
end
if nargin < 7 || isempty(outputH5FileStr)
    outputH5FileStr = [amisrRootPath,amisrFileStr(1:20),'-energyFlux.h5'];
end
if nargin < 6
    altLim = [60 200];
end

if nargin < 5
    energyBin = logspace(3,6,30);
end

if nargin < 4
    projectionAltPFISR = 60;
end

% Energy Flux from PFISR
amisrFileNameStr = [amisrRootPath,amisrFileStr];
% amisr = read_amisr(amisrFileNameStr);
% amisrData = aer_to_field_aligned_coords(amisr,projectionAltPFISR);
% amisrData = interpolate_to_field_aligned_coords(amisrData,timeMinStr,timeMaxStr);
% [dataInv, magcoords, dataInputInv] = get_2D_energy_spectra(amisrData,energyBin',...
%    timeMinStr,timeMaxStr,altLim,'magnetic');
% create_energyFlux_hdf5(dataInv,dataInputInv,outputH5FileStr,'dataInv');

% Optical data from DASC
if isempty(timeMinStr)
    timeMinStr = datestr(min(amisrData.time(1,:)));
end
if isempty(timeMaxStr)
    timeMaxStr = datestr(max(amisrData.time(1,:)));
end
create_DASC_hdf5_low_memory_v2(dascRootPath,outputH5FileStr,...
projectionAltDASC,dascMinElevation,...
timeMinStr,timeMaxStr,...
dascCalFileLat,dascCalFileLon, dascSetDownloadFlag);

% data.PFISR = dataPFISREnergyFlux;
status = ['Successfully stored the file in ',outputH5FileStr];
end

