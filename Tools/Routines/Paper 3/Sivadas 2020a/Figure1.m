%% Figure 1
% Superposed epoch time series with SML and precipitation, 
% with average phase start and stop times marked.

% Date 26 April 2021

% Initialization
if exist('prevSubstormTableFileType','var')
    if prevSubstormTableFileType ~= substormTableFileType
        clear T superMag fileNameList condition;
    end
end

tic
% Loading all relevant data variables
if ~exist('omni','var')
    omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5'; %Data from 1995
    omni = extract_omni_data(omniFile);
end

% Using Colin Forsyth's substorm phases
storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\'; 
substormList = "G:\My Drive\Research\Projects\Collaborations\Colin Forsyth\sophie_output_90.asc"; %Check sensitivity
substormTableFileType = 2; % 2 -> Colin Forsyth's data
dataChange = 1;
stormList = "G:\My Drive\Research\Projects\Collaborations\Maria Walach\Walach1981to2019_stormlist.txt";

workDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir';

% Time window around onset to process and plot data
preOnsetWindow = duration(3,0,0); 
postOnsetWindow = duration(3,0,0);
%


%%
if ~exist('prevSubstormTableFileType','var')
    prevSubstormTableFileType = 0;
end

% Constructing table of substorms while tracking PFISR experiments
if ~(exist('T','var') && prevSubstormTableFileType == substormTableFileType && dataChange) 
    % runs only when the program starts for the first time
    timeMinStr = "01 Dec 2006";
    timeMaxStr = "31 Dec 2019";
    [T, superMag] = substorm_create_table(substormList, substormTableFileType,...
      timeMinStr, timeMaxStr);
    
else
    warning('Substorm table already exists. See table T.');
end

%%
if prevSubstormTableFileType == 0 % For the first run of this routine
    
    % Cataloging stored substorms
    % THIS is a problem. PFISR data is stored based on SuperMag list. We
    % need a different method, with a list of stored PFISR files, and
    % their time arrays, and precise data must be extracted from these. 
    fileNameList = struct2cell(dir([storeDir,strcat('Processed',filesep,'*_pfisrData.h5')]));
    filePathStr = strcat(strcat(storeDir,'Processed',filesep),string(fileNameList(1,:)'));
    for i=1:1:length(filePathStr)
        storedStorm.TimeArray{:,i} = datetime(h5read(filePathStr{i},'/time'),'ConvertFrom','datenum');
        storedStorm.beg(i) =storedStorm.TimeArray{i}(1);
        storedStorm.end(i) =storedStorm.TimeArray{i}(end);
        storedStorm.filePathStr{i} = filePathStr{i};
    end
    

else  
    warning(['Since I think this is not the first run,',...
        ' did not load data variables again.',...
        ' Set prevSubstormTableFileType=0 to rerun this section.']);
end


%% Identifying and defining isolated substorms

% Identify peak SML and peak SML Time during the substorm window 
% Substorm windows - 1 hr and + 30 minutes from the onset
% (a proxy for the strength of the substorm)

if ~(exist('T','var') && prevSubstormTableFileType == substormTableFileType) 
    
    for i=1:1:length(T.stormID)
     
        if sum(storedStorm.beg<=T.Time(i) & storedStorm.end>T.Time(i))>0
            T.storageLocation(i) = filePathStr(storedStorm.beg<=T.Time(i) & storedStorm.end>T.Time(i));
        end
    
    tempTime = T.Time(i)-duration(1,0,0):duration(0,1,0):T.Time(i)+duration(0,30,0);
    [T.peakSML(i),indx2] = min(omni.Fsml(datenum(tempTime)));
    T.peakSMLTime(i)= tempTime(indx2);
    
    
    
    end

    T = identify_qualified_substorm_duration(T); 
    
else
    
    warning('Substorm table already exists. See table T.');

end

%% Loading OMNI data corresponding to each substorm onto the table

if ~(exist('T','var') && prevSubstormTableFileType == substormTableFileType) 
    T = add_omni_array_to_table(T, omni, preOnsetWindow, postOnsetWindow);
else   
    warning('Substorm table already exists, so not updating data. See table T.');
end

%% Loading PFISR data corresponding to each substorm with pfisr data

if ~(exist('T','var') && prevSubstormTableFileType == substormTableFileType) 
    T=add_pfisr_array_to_table(T);
else
    warning('Substorm table already exists, so not updating data. See table T.');
end
toc
%% Plotting

condition = T.previousSubstormDuration>duration(3,0,0)... %    
    & T.nextSubstormDuration>duration(3,0,0) ... 
    & absDiffMLT(T.MLT, T.PFISR_MLT)<2 ...
     & T.MLAT<68;
%      & ~ismissing(T.storageLocation); %    

% condition = absDiffMLT(T.MLT, T.PFISR_MLT)<2 ...
%      & T.MLAT<68;
%  
plot_superposed_epoch(T, condition, {'BzArray','EArray','smlArray','pfisrNe'});

%% Functions
prevSubstormTableFileType = substormTableFileType; 



function dMLT = absDiffMLT(a,b)
    % Source: https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles 
    normMLT = mod(a-b,24);
    dMLT = min(24-normMLT, normMLT);
end

% Creating a function that does what substorm_table.m does, i.e., outputs a
% table of substorms that are linked with PFISR and DASC data given an
% input substorm file. 

function [T, superMag] = substorm_create_table(superMagFileStr, superMagFileType,...
  timeMinStr, timeMaxStr)

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/Elements/Nithin/Data/';
    storeDir = '/media/nithin/Elements/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
else
%     error(['Not configured for this computer: ',get_computer_name]);
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
end
outputAMISRFileStr = 'amisrWebDatabase.h5';
amisrDatabaseStr = [dataDir,outputAMISRFileStr];
% dascFileStr = [storeDir,'dascDatabase.h5'];
omniFileStr = [dataDir,'omni.h5'];

% Substorms at PFISR [IMPORTANT] : Range of closeness to PFISR
% Dmlt = 2;
% Dmlat = 50; %desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat

% Poker Flat geolocation
pkrGLAT = 65.126;
pkrGLON = -147.47;
pkrh0=0.693;

% Loading database

% Create amisr database
if ~isfile(amisrDatabaseStr)
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataDir,outputAMISRFileStr]);
end
%
% Load supermag database
superMag1 = extract_data(superMagFileStr,superMagFileType);
superMag = superMag1(superMag1.phase==2,:);
superMag.stormID = (1:length(superMag.time))';

for i=1:1:length(superMag.stormID)
    if i==1
        superMag.growthPhaseStart(i) =nan;
    else
        try 
        superMag.growthPhaseStart(i)=superMag1{superMag1.time > superMag.time(i-1) & superMag1.time < superMag.time(i) & superMag1.phase==1,2};
        catch
            superMag.growthPhaseStart(i)=nan;
        end
    end
    
    if i==length(superMag.stormID)
        superMag.recoveryPhaseStart(i) = nan;
        superMag.recoveryPhaseEnd(i) = nan;
    else
        try
        superMag.recoveryPhaseStart(i)=superMag1{superMag1.time > superMag.time(i) & superMag1.time < superMag.time(i+1) & superMag1.phase==3,2};
        catch
            superMag.recoveryPhaseStart(i) = nan;
        end
        
        try
            superMag.recoveryPhaseEnd(i)=superMag1{superMag1.time > superMag.time(i) & superMag1.time < superMag.time(i+1) & superMag1.phase~=3 ,2};
        catch
            superMag.recoveryPhaseEnd(i)=nan;
        end
        
    end
        
end

% Load amisr data
amisr = extract_amisr_data(amisrDatabaseStr);

% Load omni data
omni.AE = h5read(omniFileStr,'/Indices/AE');
omni.time = unixtime2matlab(h5read(omniFileStr,'/Time'));

% Calculations
% Estimating PFISR magnetic coordinates 
[superMag.pfisrMlat,superMag.pfisrMlon,superMag.pfisrMlt] = get_magnetic_coordinates([pkrh0,pkrGLAT,pkrGLON],superMag.time(:));

% Interpolating AE index of the substorms
superMag.AE = interp1(omni.time,omni.AE,superMag.time);

%%
% Selecting substorms closest to PFISR location
% deltaMLT = absDiffMLT(superMag.pfisrMlt,superMag.mlt);
% desiredMLTIndx = abs(deltaMLT)<Dmlt;
% desiredMLATIndx = superMag.mlat<superMag.pfisrMlat+Dmlat;
% closestSubstormIndx = desiredMLTIndx & desiredMLATIndx; 

%% Adding the PFISR experiments running during the substorm time

expBCArray = barker_coded_experiments();
% Barker Coded PFISR Experiment Filter
bcFilterIndx = zeros(1,length(amisr.expId));
for iexp = 1:1:length(expBCArray)
    bcFilterIndx = bcFilterIndx|strcmp(strtrim(expBCArray(iexp)),cellstr(deblank(amisr.expName)));
end

for iStorm = 1:1:length(superMag.stormID)
    tempIndx = find_amisr_exp(superMag.time(iStorm),amisr.startTime, amisr.endTime);
    numExp(iStorm) = length(tempIndx);
    amisrIndx(iStorm)=tempIndx(1);
end
superMag.expID = repmat(string("nan"),1,length(superMag.stormID))';
superMag.expName = repmat(string("nan"),1,length(superMag.stormID))';
superMag.status = repmat(string("nan"),1,length(superMag.stormID))';
superMag.startTime = nan(1,length(superMag.stormID))';
superMag.endTime = nan(1,length(superMag.stormID))';
superMag.expBC = false(1,length(superMag.stormID))';
superMag.numberOfSimultaneousExp = nan(1,length(superMag.stormID))';

superMag.expID(~isnan(amisrIndx))=amisr.expId(amisrIndx(~isnan(amisrIndx)));
superMag.expName(~isnan(amisrIndx))=amisr.expName(amisrIndx(~isnan(amisrIndx)));
superMag.status(~isnan(amisrIndx))=amisr.status(amisrIndx(~isnan(amisrIndx)));
superMag.startTime(~isnan(amisrIndx))=amisr.startTime(amisrIndx(~isnan(amisrIndx)));
superMag.endTime(~isnan(amisrIndx))=amisr.endTime(amisrIndx(~isnan(amisrIndx)));
superMag.expBC(~isnan(amisrIndx))=bcFilterIndx(amisrIndx(~isnan(amisrIndx)));
superMag.numberOfSimultaneousExp(~isnan(amisrIndx))=numExp(amisrIndx(~isnan(amisrIndx)));


%% Create a table
T = table(superMag.datetime,superMag.stormID,...
    datetime(superMag.growthPhaseStart,'ConvertFrom','datenum'),...
    datetime(superMag.recoveryPhaseStart,'ConvertFrom','datenum'),...
    datetime(superMag.recoveryPhaseEnd,'ConvertFrom','datenum'),...
    superMag.AE, superMag.mlat,...
    superMag.mlt, superMag.pfisrMlt,...
    superMag.expID,superMag.expName,...
    datetime(superMag.startTime,'ConvertFrom','datenum'),...
    datetime(superMag.endTime,'ConvertFrom','datenum'),...
    superMag.status,...
    superMag.expBC,...
    'VariableNames',{'Time','stormID','GPStart','RPStart','RPEnd','AE','MLAT','MLT','PFISR_MLT',...
    'PFISR_ExpID','PFISR_ExpName','PFISR_startTime','PFISR_endTime',...
    'PFISR_ExpStatus','BarkerCode'});

end


% A function that extracts data from substorm/storm/SCM lists
function T1 = extract_data(loadFile, ftype)
    
    if ftype == 1
        format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
        tempData = load_ascii_files(loadFile, format, 69);
        superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        superMag.time = datenum(datestr(superMag.datetime));
        superMag.mlat = tempData{6};
        superMag.mlt = tempData{7};
    end

    if ftype == 2 %Forsyth SOPHIE substorm phases data
        format = '%4f/%2f/%2f-%2f:%2f:%2f %u %u %5.2f %5.2f';
        tempData = load_ascii_files(loadFile, format, 16);

        T.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},tempData{6});
        T.time = datenum(datestr(T.datetime));
        T.phase = tempData{7};
        T.flag = tempData{8};
        T.mlt = tempData{9};
        T.mlat = tempData{10};
    end
    
    if ftype == 3 %Walach Storm data
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f';
        tempData = load_ascii_files(loadFile, format, 6);

        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeM = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.datetimeR = datetime(tempData{1+5+5},tempData{2+5+5},tempData{3+5+5},tempData{4+5+5},tempData{5+5+5},zeros(size(tempData{5+5+5})));
        T.datetimeE = datetime(tempData{1+5+5+5},tempData{2+5+5+5},tempData{3+5+5+5},tempData{4+5+5+5},tempData{5+5+5+5},zeros(size(tempData{5+5+5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeM = datenum(datestr(T.datetimeM));
        T.timeR = datenum(datestr(T.datetimeR));
        T.timeE = datenum(datestr(T.datetimeE));
        T.symH_min = tempData{21};
    end
    
     if ftype == 4 %Walach SCM data with preceeding pre without preceeding substorms
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f';
        tempData = load_ascii_files(loadFile, format, 15);
        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeE = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
     end
    
     if ftype == 5 %Luisa, current sheet scattering
        format = ['%s %4f-%2f-%2f/%2f:%2f:%2f.%3f \t %4f-%2f-%2f/%2f:%2f:%2f.%3f',repmat(' %*s',1,44)];
        tempData = load_ascii_files(loadFile, format, 0);
        T.sc = tempData{1};
        T.datetimeI = datetime(tempData{2},tempData{3},tempData{4},tempData{5},tempData{6},tempData{7},tempData{8});
        T.datetimeE = datetime(tempData{1+8},tempData{2+8},tempData{3+8},tempData{4+8},tempData{5+8},tempData{6+8},tempData{7+8});
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
    end
    
    T1 = struct2table(T);
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
    "MSWinds23_3dt                            ";
    "MSWinds23m                               ";
    "MSWinds23hr                              ";
    "MSWinds21                                ";
    "MSWinds26.v03                          ";
    "Semeter01                              ";
    "Sporadic01                             ";
    "Sporadic02                             ";
    "Sporadic03                             ";
    "Sporadic04                             ";
    "Sporadic14                             ";
    "Sporadic15                             ";
    "Sporadic15_3dt                         ";
    "ThemisD1.v01                             "
    "MSWinds27.v01                            "];
end

% Extracting Omni data and developing gridded interpolants of the database
function omni = extract_omni_data(omniFile)
    
    %Time
    time = unixtime2matlab(h5read(omniFile,'/Time'));
    
        % Auroral electrojet indices
    SML = h5read(omniFile,'/Indices/SML');
    omni.Fsml = griddedInterpolant(time,SML);
    
    SMU = h5read(omniFile,'/Indices/SMU');
    omni.Fsmu = griddedInterpolant(time,SMU);
    
    AL = h5read(omniFile,'/Indices/AL');
    AL(AL==99999)=nan;
    omni.FAL = griddedInterpolant(time, AL);
    
    symH = h5read(omniFile,'/Indices/SYM_H');
    omni.FsymH = griddedInterpolant(time, symH);
    % Solar wind
    
        % Dynamic Pressure
    Pdyn = h5read(omniFile,'/FlowPressure');
    omni.Fp = griddedInterpolant(time,Pdyn);
    
        % IMF Bz
    BzGSM = h5read(omniFile,'/BField/BzGSM');
    omni.FBz = griddedInterpolant(time,BzGSM);
    
        % IMF By
    ByGSM = h5read(omniFile,'/BField/ByGSM');
    omni.FBy = griddedInterpolant(time,ByGSM);
    
        % IMF B magnitude
    BxGSE = h5read(omniFile,'/BField/BxGSE');
    ByGSE = h5read(omniFile,'/BField/ByGSE');
    BzGSE = h5read(omniFile,'/BField/BzGSE');
    B = (BxGSE.^2+ByGSE.^2+BzGSE.^2).^0.5;
    omni.FB = griddedInterpolant(time,B);
    
        % IMF B_T (Tangential to to GSM_x)
    BzGSM(BzGSM==0) = 0.0001;
    B_T = (ByGSM.^2 + BzGSM.^2).^0.5;
    
        % IMF Clock angle
    theta_c = wrapTo2Pi(atan2(ByGSM,BzGSM));
    theta_kl = wrapTo2Pi(atan2(B_T,BzGSM));
    omni.Ftheta = griddedInterpolant(time,theta_c);
        
        % Density
    density = h5read(omniFile,'/ProtonDensity');
    omni.Fdensity = griddedInterpolant(time,density);
        
        % Velocity
    velocity = h5read(omniFile,'/Velocity/V');
    omni.Fv = griddedInterpolant(time,velocity);
    
    
    % Solarwind - Magnetosphere Coupling
    l_0 = 7*(6371*10^3);
    
    E_kl = velocity.*B_T.*(sin(theta_kl/2)).^2; %Km nT/s
    % The “geoeffective” (or “merging”) electric field [Kan and Lee, 1979]
    omni.Fekl = griddedInterpolant(time,E_kl);
    
    E = 1e-9.*1e7.*(velocity.*10.^3).*((B.*10^-9).^2).*(sin(theta_c/2)).^4*l_0^2; %GW 
    omni.FE = griddedInterpolant(time,E);
    
    % Electric field?
    EField = h5read(omniFile,'/EField');
    omni.FEfield = griddedInterpolant(time,EField);
    
    % Magnetopause parameters from Shen et al., 1993

        % r_0 is the standoff distance in R_E
    r_0 = (11.4 + 0.013.*BzGSM).*(Pdyn).^(-1./6.6); % for Bz>=0
    r_0(BzGSM<0) = (11.4 + 0.14.*BzGSM(BzGSM<0)).*(Pdyn(BzGSM<0)).^(-1./6.6); % for Bz<0
    omni.Fr_0 = griddedInterpolant(time,r_0);
    
        % alpha is the level of tail flaring
    alpha = (0.58-0.010*BzGSM).*(1+0.01*Pdyn);
    omni.Falpha = griddedInterpolant(time,alpha);
    
    
end

function T = identify_qualified_substorm_duration(T)
    % Function to identify and define isolated substorms
    % identify_qualified_substorm_duration, 
    % first filters substorms from the database that are likely legitimate
    % removing activity that occurs in the day-side
    % secondly calculates the duration between the above qualified
    % substorms
    % Adds two columns to Table T
    %   previousSubstormDuration : The amount of time before which the
    %                              previous substorm had its onset. 
    %   nextSubstormDuration     : The amount of time after which the
    %                              next substorm will have its onset
    
    MLT_qualification = T.MLT >=16 | T.MLT<=8; 
    SML_qualification = T.peakSML <= -100; 
    
    qualified_substorms = MLT_qualification & SML_qualification; 
    T.qualifiedSubstorms = qualified_substorms;
    
    stormID = T.stormID(qualified_substorms);
    time = T.Time(qualified_substorms);
    substormDuration = repmat(duration, length(time),1);
    substormDuration(2:end) = time(2:end) - time(1:end-1);
    substormDuration(1) = nan;
    
    nextSubstormDuration = repmat(duration, length(time),1);
    nextSubstormDuration(1:end-1) = time(2:end) - time(1:end-1);
    nextSubstormDuration(end) = nan;
    
    
    T.previousSubstormDuration = repmat(duration, length(T.Time),1);
    T.previousSubstormDuration(:) = nan(length(T.Time),1);
    T.previousSubstormDuration(T.stormID(stormID)) = substormDuration;
    
    T.nextSubstormDuration = repmat(duration, length(T.Time),1);
    T.nextSubstormDuration(:) = nan(length(T.Time),1);
    T.nextSubstormDuration(T.stormID(stormID)) = nextSubstormDuration;
end



function plot_superposed_epoch(T, condition, parameters)
    
    tRange = [min(T.timeArrayRelOnset(1,:)) max(T.timeArrayRelOnset(1,:))];

    h=figure; 
    totalNo = length(parameters);
    p=create_panels(h,'totalPanelNo',totalNo,'marginbottom',10,'panelHeight',40);
    
    for i = 1:1:totalNo
        p(1,i).select();
        plot_specified_parameter(T,condition,tRange,parameters{i});
        
        if i~=totalNo
            set(gca,'XTick',{});    
        else
            xlabel('Substorm Time (t)');
        end
    end

end

function plot_epoch(T, stormID, parameters)
    
    tRange = [min(T.timeArrayRelOnset(1,:)) max(T.timeArrayRelOnset(1,:))];

    h=figure; 
    totalNo = length(parameters);
    p=create_panels(h,'totalPanelNo',totalNo,'marginbottom',10,'panelHeight',40);
    
    for i = 1:1:totalNo
        p(1,i).select();
        plot_specified_parameter_one_storm(T,stormID,tRange,parameters{i});
        
        if i~=totalNo
            set(gca,'XTick',{});    
        else
            xlabel('Substorm Time (t)');
        end
    end

end

function plot_specified_parameter(T,condition,tRange,parameterString)
    
    switch(parameterString)
        case 'BzArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'BzArray'});
            ylabel('IMF Bz [nT]')
            yRange = [-3 1];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'BzArray'}));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],yRange,[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'thetaArray'
            time_plot(T.timeArrayRelOnset(1,:),sin(T{condition,'thetaArray'}./2).^2);
            ylabel('sin^2(\Theta/2)');
            yRange = [0.4 0.8];
            ylim(yRange);
            med = nanmedian(nanmedian(sin(T{condition,'thetaArray'}./2).^2));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],[0 1],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'pArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'pArray'});
            ylabel('P_{dyn} [nPa]');
            yRange = [1.6,2.8];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'pArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[1.6 2.8],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'smlArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'smlArray'});
            ylabel('SML [nT]')
            yRange = [-400,0];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'smlArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-400 0],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
         
        case 'symHArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'symHArray'},15);
            ylabel('SYM-H [nT]')
            yRange = [-30,10];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'symHArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-30 10],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'r_0Array' 
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'r_0Array'});
            ylabel('r_0 [R_E]')
            yRange = [9.5,10.5];
            ylim(yRange);
            med = nanmedian(nanmedian(T{condition,'r_0Array'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[9 11],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');
            
        case 'alphaArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'alphaArray'});
            ylabel('\alpha [a.u.]')
            % ylim([9.5,10.5]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');

        case 'EArray'
            time_plot(T.timeArrayRelOnset(1,:),T{condition,'EArray'});
            ylabel('\epsilon(t) [GW]');
            hold on;
            med = nanmedian(nanmedian(T{condition,'BzArray'}));
            plot3(tRange,[med med] ,[1 1]);
            hold on;
            plot3([duration duration],[0 600],[1 1]);
            % set(gca,'YScale','log');
            ylim([10 400]);
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition)))],'Units','normalized');           
            
        case 'pfisrNe'
            text(0.5,0.5,1,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation))))],'Units','normalized');
           
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrNe'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^10 10^12]);
            ylabel('N_e');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        case 'pfisrEflux'
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrEflux'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrEbins'}{1});
            set(gca,'ColorScale','log','CLim',[10^8 10^12],'YScale','log');
            xlim([-3,3]);
            ylabel('\phi(E)');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        case 'pfisrHall'
            time_plot_2D(T.timeArrayRelOnset(1,:),...
                T{condition & ~ismissing(T.storageLocation),'pfisrHall'},...
                T{condition & ~ismissing(T.storageLocation),'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^-6 10^-2]);
            ylabel('\Sigma_H');
            text(0.75,1.05,['#Substorms: ',num2str(length(T.stormID(condition & ~ismissing(T.storageLocation) )))],'Units','normalized');
            
        otherwise
            error(['Invalid parameter ',parameterString]);
    end
    
end

function plot_specified_parameter_one_storm(T,stormID,tRange,parameterString)
    
    switch(parameterString)
        case 'BzArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'BzArray'});
            ylabel('IMF Bz [nT]')
            yRange = [-3 1];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'BzArray'}));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],yRange,[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'thetaArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),sin(T{stormID,'thetaArray'}./2).^2);
            ylabel('sin^2(\Theta/2)');
            yRange = [0.4 0.8];
            ylim(yRange);
            med = nanmedian(nanmedian(sin(T{stormID,'thetaArray'}./2).^2));
            hold on;
            plot3(tRange, [med med],[1 1]);
            hold on;
            plot3([duration duration],[0 1],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'pArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'pArray'});
            ylabel('P_{dyn} [nPa]');
            yRange = [1.6,2.8];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'pArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[1.6 2.8],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'smlArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'smlArray'});
            ylabel('SML [nT]')
            yRange = [-400,0];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'smlArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-400 0],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
         
        case 'symHArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'symHArray'});
            ylabel('SYM-H [nT]')
            yRange = [-30,10];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'symHArray'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[-30 10],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'r_0Array' 
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'r_0Array'});
            ylabel('r_0 [R_E]')
            yRange = [9.5,10.5];
            ylim(yRange);
            med = nanmedian(nanmedian(T{stormID,'r_0Array'}));
            hold on;
            plot3(tRange,[med med],[1 1]);
            hold on;
            plot3([duration duration],[9 11],[1 1]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'alphaArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'alphaArray'});
            ylabel('\alpha [a.u.]')
            % ylim([9.5,10.5]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');

        case 'EArray'
            time_plot_stormID(T.timeArrayRelOnset(1,:),T{stormID,'EArray'});
            ylabel('\epsilon(t) [GW]');
            hold on;
            med = nanmedian(nanmedian(T{stormID,'BzArray'}));
            plot3(tRange,[med med] ,[1 1]);
            hold on;
            plot3([duration duration],[0 600],[1 1]);
            % set(gca,'YScale','log');
            ylim([10 400]);
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            
        case 'pfisrNe'
            if ~ismissing(T.storageLocation(stormID))
                time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                    T{stormID ,'pfisrNe'}{1},...
                    T{stormID ,'pfisrAlt'}{1});
                set(gca,'ColorScale','log','CLim',[10^10 10^12]);
                ylabel('N_e');
                text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
                text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
            
        case 'pfisrEflux'
            if ~ismissing(T.storageLocation(stormID))
            time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                T{stormID ,'pfisrEflux'}{1},...
                T{stormID ,'pfisrEbins'}{1});
            set(gca,'ColorScale','log','CLim',[10^8 10^12],'YScale','log');
            xlim([-3,3]);
            ylabel('\phi(E)');
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
            
        case 'pfisrHall'
            if ~ismissing(T.storageLocation(stormID))
            time_plot_stormID_2D(T.timeArrayRelOnset(1,:),...
                T{stormID ,'pfisrHall'}{1},...
                T{stormID ,'pfisrAlt'}{1});
            set(gca,'ColorScale','log','CLim',[10^-6 10^-2]);
            ylabel('\Sigma_H');
            text(0.75,1.05,['#Substorms: ',num2str((T.stormID(stormID)))],'Units','normalized');
            text(0.10,1.05,['#Onset Time: ',datestr((T.Time(stormID)))],'Units','normalized');
            else
                text(0.5,0.5,1,['#Substorms: ',num2str((T.stormID(stormID))),' No Data'],'Units','normalized');
            end
        otherwise
            error(['Invalid parameter ',parameterString]);
    end
    
end

function time_plot_stormID(t,y)

plot3(t,y,20.*ones(1,numel(t)),'r');
view(2);

end

function time_plot(t,y, pdfBinNum)

if nargin<3
    pdfBinNum = 300;
end

% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;
my = nanmean(y);
medy = nanmedian(y);
sdev = nanstd(y);
sy = sdev./((size(y,1).^0.5));

edges = linspace(min(my)-nanmean(sdev),max(my)+nanmean(sdev),pdfBinNum);
for i=1:size(y,2) 
    N(:,i) = histcounts(y(:,i),edges,'Normalization','pdf');
end
[T,Y] = meshgrid(t,edges(1:end-1));


s=surf(T,Y,N); 
view(2); 
s.EdgeColor = 'none'; 
colormap(get_colormap('w','k'));
colorbar_thin('YLabel','PDF');
hold on;


plot3(t,my,20.*ones(1,numel(t)),'r'); hold on;
plot3(t,medy,20.*ones(1,numel(t)),'m'); hold on;
plot3(t,my+sy,20.*ones(1,numel(t)),'k'); hold on;
plot3(t,my-sy,20.*ones(1,numel(t)),'k');



end

function time_plot_stormID_2D(t,y,z)


% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;

[X,Y] = meshgrid(t,z);

zValue = y';
    
h=pcolor(hours(X),Y,zValue);
set(h,'EdgeColor','none');
shading flat;
set(gca,'Layer','top');

hold on;

%% Calculating Y axis and Color axis limits
    
	miny=min(z); maxy=max(z);
	ylim([miny maxy]);
    

    
    minz=min(min(zValue)); maxz=max(max(zValue));
    if ~isnan(minz) && ~isnan(maxz)
        caxis([minz maxz]);
    end
    colormap(get_colormap('w','k'));
    colorbar_thin();

end


function time_plot_2D(t,y,z)


% plot(t, y, 'Color', [0.9,0.9,0.9,0.5]);
% hold on;
yarr = cat(3,y{:});
my = nanmean(yarr,3);
medy = nanmedian(yarr,3);
sdev = squeeze(nanstd(permute(yarr,[3 1 2])));
sy = sdev./((size(yarr,1).^0.5));

[X,Y] = meshgrid(t,z);

zValue = my';
    
h=pcolor(hours(X),Y,zValue);
set(h,'EdgeColor','none');
shading flat;
set(gca,'Layer','top');

hold on;

%% Calculating Y axis and Color axis limits
    
	miny=min(z); maxy=max(z);
	ylim([miny maxy]);
    

    
    minz=min(min(zValue)); maxz=max(max(zValue));
    if ~isnan(minz) && ~isnan(maxz)
        caxis([minz maxz]);
    end
    colormap(get_colormap('w','k'));
    colorbar_thin();

end


function T = add_pfisr_array_to_table(T,alt)
    % Adds pfisr data, electron density, Hall, and energy flux on to the
    % table
    
    if nargin<2
        alt = 60:1:140;
    end
    storedFiles=unique(T.storageLocation(~ismissing(T.storageLocation))');
    for j=1:1:length(storedFiles)
        
        try
        pfisrNe = h5read(storedFiles(j),'/inputData/Ne');
        pfisrHall = h5read(storedFiles(j),'/conductivity/hall');
        pfisrEflux = h5read(storedFiles(j),'/energy/energyFlux');
        pfisrEbins = h5read(storedFiles(j),'/energy/energyBin'); 
        pfisrtime = h5read(storedFiles(j),'/time');
        pfisralt =  h5read(storedFiles(j),'/alt');
        
            for i = T.stormID(strcmp(T.storageLocation,storedFiles(j)))'

                pfisrRelTime = interp1(T.timeArray(i,:),T.timeArrayRelOnset(i,:),pfisrtime,'linear','extrap');

                T.pfisrNe{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrNe,alt,'linear')',T.timeArrayRelOnset(i,:));
                T.pfisrHall{i} = interp1(pfisrRelTime,interp1(pfisralt,pfisrHall,alt,'linear')',T.timeArrayRelOnset(i,:));
                T.pfisrEbins{i} = pfisrEbins;
                T.pfisrEflux{i} = interp1(pfisrRelTime,pfisrEflux',T.timeArrayRelOnset(i,:));
                T.pfisrAlt{i} = alt;

            end
    
        catch ME
             disp(['File: ',storedFiles(j)]);
             getReport(ME)
             warning('Continuing despite error');
        end
        
    end
end


function T = add_omni_array_to_table(T, omni, preOnsetWindow, postOnsetWindow)

    if nargin<4
        postOnsetWindow = duration(3,0,0);
    end
    
    if nargin<3
        preOnsetWindow = duration(3,0,0);
    end
    
    smlTime = -preOnsetWindow:duration(0,1,0):postOnsetWindow;
    T.timeArray = repmat(duration,height(T),length(smlTime));
    zeroMatrix = zeros(height(T),length(smlTime));
    T.smlArray = zeroMatrix;
    T.smuArray = zeroMatrix;
    T.ALArray = zeroMatrix;
    T.BzArray = zeroMatrix;
    T.pArray = zeroMatrix;
    T.densityArray = zeroMatrix;
    T.vArray = zeroMatrix;
    T.EklArray = zeroMatrix;
    T.EArray = zeroMatrix;
    T.thetaArray = zeroMatrix;
    T.r_0Array = zeroMatrix;
    T.alphaArray = zeroMatrix;
    T.symHArray = zeroMatrix;
    
    for i=1:1:height(T)
    
        TMatrix(i,:) = datenum(T.Time(i)-preOnsetWindow : duration(0,1,0) : ...
        T.Time(i)+postOnsetWindow);
    
    end
    
    T.timeArrayRelOnset = repmat(smlTime,height(T),1);
    T.timeArray = TMatrix;
    
    T.smlArray      = omni.Fsml(TMatrix);
    T.smuArray      = omni.Fsmu(TMatrix);
    T.pArray        = omni.Fp(TMatrix); 
    T.BzArray       = omni.FBz(TMatrix); 
    T.ALArray       = omni.FAL(TMatrix); 
    T.densityArray  = omni.Fdensity(TMatrix); 
    T.vArray        = omni.Fv(TMatrix); 
    T.eklArray      = omni.Fekl(TMatrix); 
    T.EArray        = omni.FE(TMatrix); 
    T.thetaArray    = omni.Ftheta(TMatrix); 
    T.r_0Array      = omni.Fr_0(TMatrix); 
    T.alphaArray    = omni.Falpha(TMatrix);
    T.symHArray     = omni.FsymH(TMatrix);

end