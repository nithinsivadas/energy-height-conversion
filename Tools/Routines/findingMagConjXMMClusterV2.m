%%
[rootPathStr, dataPathStr] = initialize_root_path;
addpath([rootPathStr,'Geopack2008']);
if isunix, setenv('LD_LIBRARY_PATH',''); end
% magFieldModel=7; %TS96 
% xmmOutputFileStr = 'xmm_eq_magnetic_point_TS96 - Copy.dat';
% c4OutputFileStr = 'c4_eq_magnetic_point_TS96 - Copy.dat';
% parobj=parpool;
% 
magFieldModel=9; %TS01 
xmmOutputFileStr = 'xmm_eq_magnetic_point_TS01.dat';
c4OutputFileStr = 'c4_eq_magnetic_point_TS01.dat';

if ispc
    xmmDataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
    c4DataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
    omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
elseif isunix
    xmmDataFolder = [initialize_root_path,'LargeFiles/issiteam2017/'];
    c4DataFolder = [initialize_root_path,'LargeFiles/issiteam2017/'];
    omniH5FileStr = [initialize_root_path,'LargeFiles/omni/omni.h5'];
end

xmmOutputFileStr = [dataPathStr,xmmOutputFileStr];
c4OutputFileStr = [dataPathStr,c4OutputFileStr];

%% Inputs
% Cropping time
timeMinXMM = datenum('2000-02-16 21:16:00');
timeMaxXMM = datenum('2012-08-30 09:29:00');
timeMinC4 = datenum('2001-01-09 15:23:00');
timeMaxC4 = datenum('2015-04-30 23:59:00');
% 
%timeMinXMM = datenum('2001-01-09 15:23:00');
%timeMaxXMM = datenum('2001-01-10 09:29:00');
%timeMinC4 = datenum('2001-01-09 15:23:00');
%timeMaxC4 = datenum('2001-01-10 09:29:00');


%%
c4MatFile = [c4DataFolder,'Mat',filesep,'cluster_all_records.mat'];
if ~isfile(c4MatFile)
    c4CoordFile = 'cluster_omni_for_t96_allrecords.txt';
    if ~isfile([c4DataFolder,c4CoordFile,'.bak'])
        [s1,msg1] = replaceinfile(',','',[c4DataFolder,c4CoordFile]);
    end
    clusterA = importdata([c4DataFolder,c4CoordFile]);
    c4.time = datenum(clusterA.textdata(2:end,1),'yyyy-mm-dd HH:MM:ss');
    c4.xGSE  = clusterA.data(:,1:3);
    c4.AE = clusterA.data(:,4);
    c4.Dst = clusterA.data(:,5);
    c4.BYimf = clusterA.data(:,6);
    c4.BZimf = clusterA.data(:,7);
    c4.Pdyn = clusterA.data(:,8);
    save(c4MatFile,'c4*');
else
    load(c4MatFile);
end
disp('1/7 Done loading cluster4 coordinates...');
%%


xmmMatFile = [xmmDataFolder,'Mat',filesep,'xmm_all_records.mat'];
if ~isfile(xmmMatFile)
    xmmCoordFile = 'xmm_omni_for_t96_allrecords.txt';
    if ~isfile([xmmDataFolder,xmmCoordFile,'.bak'])
    [s2,msg2] = replaceinfile(',','',[xmmDataFolder,xmmCoordFile]);
    end
    XMMA = importdata([xmmDataFolder,xmmCoordFile]);
    xmm.time = datenum(XMMA.textdata(2:end,1),'yyyy-mm-dd HH:MM:ss');
    xmm.xGSE  = XMMA.data(:,1:3);
    xmm.AE = XMMA.data(:,4);
    xmm.Dst = XMMA.data(:,5);
    xmm.BYimf = XMMA.data(:,6);
    xmm.BZimf = XMMA.data(:,7);
    xmm.Pdyn = XMMA.data(:,8);
    save(xmmMatFile,'xmm*');
else
    load(xmmMatFile);
end
disp('2/7 Done loading XMM coordinates...');
%%
[spacecraft.c4.xGSE,spacecraft.c4.time] = crop_time(c4.xGSE,c4.time,timeMinC4,timeMaxC4);
[spacecraft.xmm.xGSE,spacecraft.xmm.time] = crop_time(xmm.xGSE,xmm.time,timeMinXMM,timeMaxXMM);
time = spacecraft.xmm.time;

% Conversion from GSE to GEO
c=define_universal_constants;
% RE = c.RE/1000; %in km
spacecraft.c4.GEO   = onera_desp_lib_rotate(spacecraft.c4.xGSE,'gse2geo',spacecraft.c4.time);
spacecraft.xmm.GEO   = onera_desp_lib_rotate(spacecraft.xmm.xGSE,'gse2geo',spacecraft.xmm.time);

c4.maginput = zeros(length(c4.time),25);
% c4.maginput(:,2) =c4.Dst; 
% c4.maginput(:,5) = c4.Pdyn;
% c4.maginput(:,6) = c4.BYimf; 
% c4.maginput(:,7) = c4.BZimf;
[spacecraft.c4.maginput,spacecraft.c4.time] = crop_time(c4.maginput,c4.time,timeMinC4,timeMaxC4);
% xmm.Pdyn(1:5) = [1.93; 1.16; 1.15; 1.35; nan];
xmm.maginput = zeros(length(xmm.time),25);
% xmm.maginput(:,2) =xmm.Dst; 
% xmm.maginput(:,5) = xmm.Pdyn;
% BimfXMM = onera_desp_lib_rotate([zeros(length(xmm.BYimf),1), xmm.BYimf, xmm.BZimf], 'gse2gsm', xmm.time);
% xmm.maginput(:,6) = BimfXMM(:,2); 
% xmm.maginput(:,7) = BimfXMM(:,3);
[spacecraft.xmm.maginput,spacecraft.xmm.time] = crop_time(xmm.maginput,xmm.time,timeMinXMM,timeMaxXMM);
disp('3/7 Done converting coordinates from GSE to GEO...');
%% OMNI Maginput
[omnimaginput,omniTime] = generate_maginput(omniH5FileStr, timeMinXMM, timeMaxC4);
% nMonth = months(datestr(timeMin,'mmm dd yyyy'),datestr(timeMax,'mmm dd yyyy'))+1;
% dateVec = datevec(timeMin);
% dateVec(3) = 1; %setting day to 1st
% omniDataProcessTimes = datenum(dateVec);
% omniDataStorage='G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
% omnimaginput = [];
% omniTime = [];
% multiWaitbar('Generate OMNI Input',0);
% dmonth = 1./nMonth;
%     for iMonth = 1:1:nMonth
%         multiWaitbar('Generate OMNI Input','Increment',dmonth);
%         omniDataTemp = process_omni_data(datestr(omniDataProcessTimes,'yyyy mm dd'),omniDataStorage,true);
%         omnimaginput = [omnimaginput; omniDataTemp.minutely.maginput];
%         omniTime = [omniTime; omniDataTemp.minutely.time'];
%         omniDataProcessTimes = addtodate(omniDataProcessTimes,1,'month');
%     end
% multiWaitbar('Generate OMNI Input',1);
% omniPdyn = omnimaginput(:,5);
% omniPdyn(omniPdyn>99.9) = nan;
% omnimaginput(:,5) = interp_nans(omniPdyn);
omnimaginput = filter_irbem_maginput(magFieldModel,omnimaginput);
omnimaginput = interp_nans(omnimaginput);
xmmMaginput = interp1(omniTime,omnimaginput,spacecraft.xmm.time);
c4Maginput = interp1(omniTime,omnimaginput,spacecraft.c4.time);
spacecraft.xmm.maginput = xmmMaginput;
spacecraft.c4.maginput = c4Maginput;
disp('4/7 Done downloading OMNI data...');
%%
options = [0,0,0,0,0];
sysaxes = 1; %GEO Input coordinates
stopAlt = 110;

%%
% multiWaitbar('XMM Calculation');
% dt = 1./length(spacecraft.xmm.time);
tic
for i = 1:1:length(spacecraft.xmm.time)
    
    try
    [foot, BGEO, eqGEO, eqBGEO] = nithin_find_foot_point_fast(magFieldModel,100,...
    sysaxes,spacecraft.xmm.time(i),spacecraft.xmm.GEO(i,1),...
    spacecraft.xmm.GEO(i,2),spacecraft.xmm.GEO(i,3),...
    stopAlt,spacecraft.xmm.maginput(i,:));
    spacecraft.xmm.GDZNorth(i,:) = foot.north;
    spacecraft.xmm.GDZSouth(i,:) = foot.south;
    spacecraft.xmm.footType(i) = foot.number;
    spacecraft.xmm.eqGEO(i,:) = eqGEO;
    spacecraft.xmm.eqBGEO(i,:) = eqBGEO;
    spacecraft.xmm.BGEO(i,:) = BGEO;
    
    catch
        spacecraft.xmm.GDZNorth(i,:) = nan(1,3);
        spacecraft.xmm.GDZSouth(i,:) = nan(1,3);
        spacecraft.xmm.footType(i) = 3;
        spacecraft.xmm.eqGEO(i,:) = nan(1,3);
        spacecraft.xmm.eqBGEO(i,:) = nan(1,3);
    end
%     multiWaitbar('XMM Calculation','Increment',dt);
end
toc
spacecraft.xmm.BGSE = onera_desp_lib_rotate(spacecraft.xmm.BGEO,'geo2gse',spacecraft.xmm.time);
spacecraft.xmm.eqGSE = onera_desp_lib_rotate(spacecraft.xmm.eqGEO,'geo2gse',spacecraft.xmm.time);
spacecraft.xmm.eqBGSE = onera_desp_lib_rotate(spacecraft.xmm.eqBGEO,'geo2gse',spacecraft.xmm.time);

fileID = fopen(xmmOutputFileStr,'w');
fprintf(fileID,'%s\n','DateTime FootType Bx_GSE By_GSE Bz_GSE Eq_x_GSE Eq_y_GSE Eq_z_GSE Beq_x_GSE Beq_y_GSE Beq_z_GSE');
kipsFormat = "%19s\t%1u\t%6.2f\t%6.2f\t%6.2f\t%7.2f\t%7.2f\t%7.2f\t%6.2f\t%6.2f\t%6.2f\n";

for i=1:1:length(spacecraft.xmm.time)
    fprintf(fileID,kipsFormat,...
        (datestr(spacecraft.xmm.time(i),'yyyy-mm-dd HH:MM:ss'))',...
        spacecraft.xmm.footType(i),...    
        spacecraft.xmm.BGSE(i,:),...
        spacecraft.xmm.eqGSE(i,:),...
        spacecraft.xmm.eqBGSE(i,:));
%         spacecraft.xmm.maginput(i,1:9),...
%         spacecraft.xmm.xGSE(i,:),...
    %         spacecraft.xmm.GDZNorth(i,:),...
%         spacecraft.xmm.GDZSouth(i,:));
        
end    
fclose(fileID);

disp('Done writing XMM magnetic coordinates output data file...');

disp('5/7 Done calculating XMM...');
%% Cluster Calculation
%multiWaitbar('Cluster Calculation');
% dt = 1./length(spacecraft.c4.time);
tic
for i = 1:1:length(spacecraft.c4.time)
    
    try
    [foot, BGEO, eqGEO, eqBGEO] = nithin_find_foot_point_fast(magFieldModel,100,...
    sysaxes,spacecraft.c4.time(i),spacecraft.c4.GEO(i,1),...
    spacecraft.c4.GEO(i,2),spacecraft.c4.GEO(i,3),...
    stopAlt,spacecraft.c4.maginput(i,:));
    spacecraft.c4.GDZNorth(i,:) = foot.north;
    spacecraft.c4.GDZSouth(i,:) = foot.south;
    spacecraft.c4.footType(i) = foot.number;
    spacecraft.c4.eqGEO(i,:) = eqGEO;
    spacecraft.c4.eqBGEO(i,:) = eqBGEO;
    spacecraft.c4.BGEO(i,:) = BGEO;
    
    catch
        spacecraft.c4.GDZNorth(i,:) = nan(1,3);
        spacecraft.c4.GDZSouth(i,:) = nan(1,3);
        spacecraft.c4.footType(i) = 3;
        spacecraft.c4.eqGEO(i,:) = nan(1,3);
        spacecraft.c4.eqBGEO(i,:) = nan(1,3);
    end
%    multiWaitbar('Cluster Calculation','Increment',dt);
end
% multiWaitbar('CloseAll');
toc
spacecraft.c4.BGSE = onera_desp_lib_rotate(spacecraft.c4.BGEO,'geo2gse',spacecraft.c4.time);
spacecraft.c4.eqGSE = onera_desp_lib_rotate(spacecraft.c4.eqGEO,'geo2gse',spacecraft.c4.time);
spacecraft.c4.eqBGSE = onera_desp_lib_rotate(spacecraft.c4.eqBGEO,'geo2gse',spacecraft.c4.time);


% Writing to an ASCII File


% celltable =[cellstr(datestr(spacecraft.xmm.time,'yyyy-mm-dd HH:MM:ss')),...
%     num2cell(spacecraft.xmm.xGSE(:,1)),num2cell(spacecraft.xmm.xGSE(:,2)),num2cell(spacecraft.xmm.xGSE(:,3)),...
%     num2cell(round(spacecraft.xmm.BGSE(:,1),2)),num2cell(round(spacecraft.xmm.BGSE(:,2),2)),num2cell(round(spacecraft.xmm.BGSE(:,3),2)),...
%     num2cell(spacecraft.xmm.footType),...
%     num2cell(round(spacecraft.xmm.eqGSE(:,1),2)),num2cell(round(spacecraft.xmm.eqGSE(:,2),2)),num2cell(round(spacecraft.xmm.eqGSE(:,3),2)),...
%     num2cell(round(spacecraft.xmm.GDZNorth(:,1),2)),num2cell(round(spacecraft.xmm.GDZNorth(:,2),2)),num2cell(round(spacecraft.xmm.GDZNorth(:,3),2)),...
%     num2cell(round(spacecraft.xmm.GDZSouth(:,1),2)),num2cell(round(spacecraft.xmm.GDZSouth(:,2),2)),num2cell(round(spacecraft.xmm.GDZSouth(:,3),2))];
% T = cell2table(celltable,'VariableNames',{'DateTime', 'x_GSE_RE', 'y_GSE_RE', 'z_GSE_RE',...
%     'Bx_GSE_nT','By_GSE_nT','Bz_GSE_nT','footType','x_eqGSE_RE','y_eqGSE_RE','z_eqGSE_RE',...
%     'alt_NFoot_km','lat_NFoot_deg','East_lon_NFoot_deg',...
%     'alt_SFoot_km','lat_SFoot_deg','East_lon_SFoot_deg'});
% writetable(T,'tabledata.dat','Delimiter','\t');

% formatSpec = "%s %8u %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8u %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n";
% fprintf(fileID,'%15s %12s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n','DateTime',...
%     'Kp','Dst','Np_SW','V_SW','Pdyn','yIMF_GSM','zIMF_GSM','G1','G2',...
%     'x_GSE', 'y_GSE', 'z_GSE',...
%     'Bx_GSE','By_GSE','Bz_GSE','footType','x_eqGSE','y_eqGSE','z_eqGSE',...
%     'alt_NFoot','lat_NFoot','lon_NFoot',...
%     'alt_SFoot','lat_SFoot','lon_SFoot');
% fprintf(fileID,'%5s %7s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n','[yyyy-mm-dd HH:MM:ss]',...
%     '[a.u]','[nT]','[cm-3]','[km/s]','[nPa]','[nT]','[nT]','[a.u.]','[a.u.]',...
%     '[RE]', '[RE]', '[RE]',...
%     '[nT]','[nT]','[nT]','[0,1,2]','[RE]','[RE]','[RE]',...
%     '[km]','[deg]','[deg]',...
%     '[km]','[deg]','[deg]');

fileID = fopen(c4OutputFileStr,'w');
% formatSpec = "%s %8u %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8u %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n";
% fprintf(fileID,'%15s %12s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n','DateTime',...
%     'Kp','Dst','Np_SW','V_SW','Pdyn','yIMF_GSM','zIMF_GSM','G1','G2',...
%     'x_GSE', 'y_GSE', 'z_GSE',...
%     'Bx_GSE','By_GSE','Bz_GSE','footType','x_eqGSE','y_eqGSE','z_eqGSE',...
%     'alt_NFoot','lat_NFoot','lon_NFoot',...
%     'alt_SFoot','lat_SFoot','lon_SFoot');
% fprintf(fileID,'%5s %7s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n','[yyyy-mm-dd HH:MM:ss]',...
%     '[a.u]','[nT]','[cm-3]','[km/s]','[nPa]','[nT]','[nT]','[a.u.]','[a.u.]',...
%     '[RE]', '[RE]', '[RE]',...
%     '[nT]','[nT]','[nT]','[0,1,2]','[RE]','[RE]','[RE]',...
%     '[km]','[deg]','[deg]',...
%     '[km]','[deg]','[deg]');
fprintf(fileID,'%s \n','DateTime FootType Bx_GSE By_GSE Bz_GSE Eq_x_GSE Eq_y_GSE Eq_z_GSE Beq_x_GSE Beq_y_GSE Beq_z_GSE');

for i=1:1:length(spacecraft.c4.time)
        fprintf(fileID,kipsFormat,...
        (datestr(spacecraft.c4.time(i),'yyyy-mm-dd HH:MM:ss'))',...
        spacecraft.c4.footType(i),...    
        spacecraft.c4.BGSE(i,:),...
        spacecraft.c4.eqGSE(i,:),...
        spacecraft.c4.eqBGSE(i,:));
    
%     fprintf(fileID,formatSpec,...
%         (datestr(spacecraft.c4.time(i),'yyyy-mm-dd HH:MM:ss'))',...
%         spacecraft.c4.maginput(i,1:9),...
%         spacecraft.c4.xGSE(i,:),...
%         spacecraft.c4.BGSE(i,:),...
%         spacecraft.c4.footType(i),...
%         spacecraft.c4.eqGSE(i,:),...
%         spacecraft.c4.GDZNorth(i,:),...
%         spacecraft.c4.GDZSouth(i,:));
end    
fclose(fileID);
disp('Done writing C4 magnetic coordinates output data file...');

disp('6/7 Done calculating Cluster...');
disp('7/7 Process Complete');
%%
% conjunction.probeStr='xmm';
% conjunction.radius = 500;
% 
% probes = find_magnetic_conjunctions_at_mag_equator(time,spacecraft,conjunction,110,1,magFieldModel);
% probes.time = time;
% multiWaitbar('closeall');

%% Generate maginput from omni.h5
function [maginput,time] = generate_maginput(omniH5FileStr, timeMin, timeMax)
    omniTime = unixtime2matlab(h5read(omniH5FileStr,'/Time'));
    omniIndx = 1:1:length(omniTime);
    minTimeIndx = find_time(omniTime,datestr(timeMin));
    maxTimeIndx = find_time(omniTime,datestr(timeMax));
    deltaTimeIndx = maxTimeIndx - minTimeIndx +1;
    timeIndx = minTimeIndx:1:maxTimeIndx;
    
    maginput(:,1) = h5read(omniH5FileStr, '/Indices/Kp', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,2) = h5read(omniH5FileStr, '/Indices/SYM_H', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,3) = h5read(omniH5FileStr, '/ProtonDensity', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,4) = h5read(omniH5FileStr, '/Velocity/V', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,5) = h5read(omniH5FileStr, '/FlowPressure', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,6) = h5read(omniH5FileStr, '/BField/ByGSM', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,7) = h5read(omniH5FileStr, '/BField/BzGSM', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,8) = h5read(omniH5FileStr, '/TSY/G1', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,9) = h5read(omniH5FileStr, '/TSY/G2', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,10) = h5read(omniH5FileStr, '/TSY/G3', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,11) = h5read(omniH5FileStr, '/TSY/W1', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,12) = h5read(omniH5FileStr, '/TSY/W2', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,13) = h5read(omniH5FileStr, '/TSY/W3', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,14) = h5read(omniH5FileStr, '/TSY/W4', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,15) = h5read(omniH5FileStr, '/TSY/W5', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,16) = h5read(omniH5FileStr, '/TSY/W6', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,17) = h5read(omniH5FileStr, '/Indices/AL', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,18:25) = nan(length(timeIndx),8);
    
    time = omniTime(timeIndx);
end

%% Replace string in file
function [s, msg] = replaceinfile(str1, str2, infile, outfile)
%REPLACEINFILE replaces characters in ASCII file using PERL
%  
% [s, msg] = replaceinfile(str1, str2, infile)
%    replaces str1 with str2 in infile, original file is saved as "infile.bak"
%
% [s, msg] = replaceinfile(str1, str2, infile, outfile)
%    writes contents of infile to outfile, str1 replaced with str2
%    NOTE! if outputfile is '-nobak' the backup file will be deleted
%
% [s, msg] = replaceinfile(str1, str2)
%    opens gui for the infile, replaces str1 with str2 in infile, original file is saved as "infile.bak"
%
% in:  str1      string to be replaced
%      str2      string to replace with
%      infile    file to search in
%      outfile   outputfile (optional) if '-nobak'
%
% out: s         status information, 0 if succesful
%      msg       messages from calling PERL 

% Pekka Kumpulainen 30.08.2000
% 16.11.2008 fixed for paths having whitespaces, 
% 16.11.2008 dos rename replaced by "movefile" to force overwrite
% 08.01.2009 '-nobak' option to remove backup file, fixed help a little..
%
% TAMPERE UNIVERSITY OF TECHNOLOGY  
% Measurement and Information Technology
% www.mit.tut.fi

message = nargchk(2,4,nargin);
if ~isempty(message)
    error(message)
end

%% check inputs
if ~(ischar(str1) && ischar(str2))
    error('Invalid string arguments.')
end
% in case of single characters, escape special characters 
% (at least someof them)
switch str1
    case {'\' '.'}
        str1 = ['\' str1];
end

%% uigetfile if none given
if nargin < 3;
    [fn, fpath] = uigetfile('*.*','Select file');
    if ~ischar(fn)
        return
    end
    infile = fullfile(fpath,fn);
end

%% The PERL stuff
perlCmd = sprintf('"%s"',fullfile(matlabroot, 'sys\perl\win32\bin\perl'));
perlstr = sprintf('%s -i.bak -pe"s/%s/%s/g" "%s"', perlCmd, str1, str2,infile);

[s,msg] = dos(perlstr);

%% rename files if outputfile given
if ~isempty(msg)
    error(msg)
else
    if nargin > 3 % rename files
        if strcmp('-nobak',outfile)
            delete(sprintf('%s.bak',infile));
        else
            movefile(infile, outfile);
            movefile(sprintf('%s.bak',infile), infile);
        end
    end
end
end
