%% Substorms in the vicinity of PFISR
% Find SuperMag Substorms and PFISR conjunctions
clear all;

%% Initialization
dataDir = 'G:\My Drive\Research\Projects\Data\';
outputFileStr = 'amisrWebDatabase.h5';
amisrDatabaseStr = [dataDir,outputFileStr];

superMagFileStr = 'G:\My Drive\Research\Projects\Paper 3\Data\substorms_superMag_20190530.txt';

omniFileStr = [dataDir,'omni.h5'];

timeMinStr = "01 Dec 2006";
timeMaxStr = "31 Jul 2019";

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
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataDir,outputFileStr]);
end
%%
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
%% Substorms AE>500 nT, 2<MLT<6

filterIndx = superMag.mlt<8 & superMag.mlt>5 & superMag.AE>500 & superMag.AE<3000;

%%
% figure; 
% plot_polar_scatter(superMag.mlat(closestSubstormIndx),superMag.mlt(closestSubstormIndx),'RLim',[40,90]);
% %%
% figure; 
% plot_polar_scatter(superMag.mlat,superMag.mlt,'RLim',[40,90]);
% 
% %%
% figure; 
% plot_polar_scatter(superMag.mlat,superMag.mlt,'RLim',[40,90],'zColor',superMag.AE,'zLim',[0,1200],'markerSize',2);
% %%
% figure; 
% plot_polar_scatter(superMag.mlat(closestSubstormIndx),superMag.mlt(closestSubstormIndx),'RLim',[60,70],'zColor',superMag.AE(closestSubstormIndx),'zLim',[0,600],'markerSize',2);

%% Adding the PFISR experiments running during the substorm time

expBCArray = barker_coded_experiments();
% Filter of Barker Coded PFISR Experiments
bcFilterIndx = zeros(1,length(amisr.expId));
for iexp = 1:1:length(expBCArray)
    bcFilterIndx = bcFilterIndx|strcmp(strtrim(expBCArray(iexp)),cellstr(deblank(amisr.expName)));
end

for iStorm = 1:1:length(superMag.stormID)
    amisrIndx = find_amisr_exp(superMag.time(iStorm),amisr.startTime, amisr.endTime);
    if ~isnan(amisrIndx)
        superMag.expID(iStorm) = amisr.expId(amisrIndx(1));
        superMag.expName(iStorm) = amisr.expName(amisrIndx(1));
        superMag.status(iStorm) = amisr.status(amisrIndx(1));
        superMag.startTime(iStorm) = amisr.startTime(amisrIndx(1));
        superMag.endTimeTime(iStorm) = amisr.endTime(amisrIndx(1));
        superMag.additionalComments(iStorm) = {' '};
        superMag.expBC(iStorm) = bcFilterIndx(amisrIndx(1));
              
        
        if length(amisrIndx)>1
            superMag.additionalComments(iStorm,:) = {['More than ',num2str(length(amisrIndx)),' experiments were running']};
        end 
    else
        superMag.expID(iStorm) = "nan";
        superMag.expName(iStorm) = "nan";
        superMag.status(iStorm) = "nan";
        superMag.startTime(iStorm) = nan;
        superMag.endTimeTime(iStorm) = nan;
        superMag.expBC(iStorm) = false;
    end
end

%% Create a table
T = table(superMag.datetime(closestSubstormIndx),...
    superMag.AE(closestSubstormIndx), superMag.mlat(closestSubstormIndx),...
    superMag.mlt(closestSubstormIndx),...
    superMag.expID(closestSubstormIndx)',superMag.expName(closestSubstormIndx)',...
    superMag.status(closestSubstormIndx)',...
    superMag.expBC(closestSubstormIndx)',...
    'VariableNames',{'Time','AE','MLAT','MLT','PFISR_ExpID','PFISR_ExpName','PFISR_ExpStatus','BarkerCode'});
%%
disp('Table of Barker Code Experiments during a SuperMag substorm');
T(T.BarkerCode,:)
%%
disp('Table of PFISR Experiments during a SuperMag substorm');
T(~strcmp(T.PFISR_ExpID,"nan"),:)

%% Finding if DASC is ON during a particular substorm with PFISR ON
substormIndx = 1:1:length(superMag.time);
cIndx = substormIndx(closestSubstormIndx);
 for cStorm = 1:1:length(cIndx)
        superMagCloseStorm.time(cStorm) = superMag.time(cIndx(cStorm));
        [timeStamp, wavelength] = find_DASC_FITS(datestr(superMagCloseStorm.time(cStorm)));
        if ~isnan(timeStamp)
            superMagCloseStorm.DASC_timeMin(cStorm) = min(timeStamp);
            superMagCloseStorm.DASC_timeMax(cStorm) = max(timeStamp);
            superMagCloseStorm.DASC_wavelength(cStorm) = {num2str(unique(wavelength))};
        else
            superMagCloseStorm.DASC_timeMin(cStorm) = nan;
            superMagCloseStorm.DASC_timeMax(cStorm) = nan;
            superMagCloseStorm.DASC_wavelength(cStorm) = {'nan'};
        end
 end
 %%
 for cStorm = 1:1:length(cIndx)
     if ~isnan(superMagCloseStorm.DASC_timeMin)
        superMagCloseStorm.DASC_timeMin1(cStorm) = datetime(superMagCloseStorm.DASC_timeMin(cStorm));
     end
 end
 
 %%
 T1 = table(superMag.datetime(closestSubstormIndx),...
    superMag.AE(closestSubstormIndx), superMag.mlat(closestSubstormIndx),...
    superMag.mlt(closestSubstormIndx),...
    superMag.expID(closestSubstormIndx)',superMag.expName(closestSubstormIndx)',...
    superMag.status(closestSubstormIndx)',...
    superMag.expBC(closestSubstormIndx)',...
    datetime(superMagCloseStorm.DASC_timeMin','ConvertFrom','datenum'),...
    datetime(superMagCloseStorm.DASC_timeMax','ConvertFrom','datenum'),...
    superMagCloseStorm.DASC_wavelength',...
    'VariableNames',...
    {'Time','AE','MLAT','MLT','PFISR_ExpID',...
    'PFISR_ExpName','PFISR_ExpStatus','BarkerCode',...
    'DASC_TimeMin','DASC_TimeMax','DASC_Wavelength'});
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

function [amisrIndx] = find_amisr_exp(time, startTimeArr, endTimeArr)
% Finds the amisr experiment indx
    amisrIndx = find(time>startTimeArr & time<endTimeArr);
    if isempty(amisrIndx)
        amisrIndx=nan;
    end
    
end

function expArr = barker_coded_experiments()
expArr =["GenPINOT_PulsatingAurora_TN30          ";
    "Inspire_v01                            ";
    "Kelley01                               ";
    "MSWinds23                              ";
    "MSWinds23_dt013                        ";
    "MSWinds26.v03                          ";
    "Semeter01                              ";
    "Sporadic04                             ";
    "ThemisD1.v01                             "];
 end