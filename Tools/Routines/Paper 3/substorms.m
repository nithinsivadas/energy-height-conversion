%% Substorms in the vicinity of PFISR
% Find SuperMag Substorms and PFISR conjunctions
clear all;

%% Initialization
dataStoreDir = "G:\Team Drives\Semeter-Research in Progress\All AMISR Experiments\";
outputFileStr = "amisrWebDatabase_temp.h5";
amisrDatabaseStr = strcat(dataStoreDir,outputFileStr);

superMagFileStr = "G:\My Drive\Research\Projects\Paper 3\Data\substorms_superMag_20190530.txt";

omniFileStr = "G:\My Drive\Research\Projects\Data\omni.h5";

timeMinStr = "01 May 2018";
timeMaxStr = "01 Dec 2018";

% Substorms at PFISR [IMPORTANT]
Dmlt = 2;
Dmlat = 1;

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

%% Loading database
% Create amisr database
if ~isfile(amisrDatabaseStr)
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataStoreDir,outputFileStr]);
end

% Load supermag database
superMag = extract_superMag_data(superMagFileStr);
superMag.stormID = (1:length(superMag.time))';

% Load amisr data
amisr = extract_amisr_data(amisrDatabaseStr);
%% Load omni data
omni.AE = h5read(omniFileStr,'/Indices/AE');
omni.time = unixtime2matlab(h5read(omniFileStr,'/Time'));
%% Calculations
% Estimating PFISR magnetic coordinates 
[superMag.pfisrMlat,superMag.pfisrMlon,superMag.pfisrMlt] = get_magnetic_coordinates([pkrh0,pkrGLAT,pkrGLON],superMag.time(:));
%% Selecting AE index of the substorms
superMag.AE = interp1(omni.time,omni.AE,superMag.time);

%%
% Selecting substorms closest to PFISR location
deltaMLT = mod(superMag.pfisrMlt - superMag.mlt,24);
desiredMLTIndx = abs(deltaMLT)<Dmlt;
desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat;
closestSubstormIndx = desiredMLTIndx & desiredMLATIndx; 


%%
figure; 
plot_polar_scatter(superMag.mlat(closestSubstormIndx),superMag.mlt(closestSubstormIndx),'RLim',[40,90]);
%%
figure; 
plot_polar_scatter(superMag.mlat,superMag.mlt,'RLim',[40,90]);

%%
figure; 
plot_polar_scatter(superMag.mlat,superMag.mlt,'RLim',[40,90],'zColor',superMag.AE,'zLim',[0,600],'markerSize',2);
%%
figure; 
plot_polar_scatter(superMag.mlat(closestSubstormIndx),superMag.mlt(closestSubstormIndx),'RLim',[60,70],'zColor',superMag.AE(closestSubstormIndx),'zLim',[0,600],'markerSize',2);
%% Functions
function data=load_ascii_files(loadFile, format, headerlines)
if nargin<3
    headerlines = 1;
end
fileID = fopen(loadFile,'r');
data = textscan(fileID, format, 'headerlines', headerlines);
fclose(fileID);
end

function superMag = extract_superMag_data(loadFile)
format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
tempData = load_ascii_files(loadFile, format, 69);
superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
superMag.time = datenum(datestr(superMag.datetime));
superMag.mlat = tempData{6};
superMag.mlt = tempData{7};
end

function amisr = extract_amisr_data(loadFile)
    table = read_h5_data(loadFile);
    amisr.expId = string(table.Data{2});
    amisr.expName = string(table.Data{3});
    amisr.status = string(table.Data{4});
    amisr.startTime = unixtime2matlab(table.Data{5});
    amisr.endTime = unixtime2matlab(table.Data{1});
end

function [Lm,MLT] = get_pfisr_magnetic_coordinates(time,maginput,GDZ,magFieldNo)
    [Lm,~,~,~,~,MLT] = onera_desp_lib_make_lstar(magFieldNo,[0,0,0,0,0],0,time,GDZ(3),GDZ(1),GDZ(2),maginput);
%     tup=py.aacgmv2.wrapper.get_aacgm_coord(GDZ(1), GDZ(2), GDZ(3), time, 'TRACE');
%     MLAT = double(py.array.array('d',py.numpy.nditer(tup{1})));
%     MLON = double(py.array.array('d',py.numpy.nditer(tup{2})));
%     MLT_AACGM = double(py.array.array('d',py.numpy.nditer(tup{3})));
    Lm = abs(Lm);
end