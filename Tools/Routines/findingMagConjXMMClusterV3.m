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
xmmOutputFileStr = 'xmm_eq_magnetic_point_TS01-copy.dat';
c4OutputFileStr = 'c4_eq_magnetic_point_TS01-copy.dat';

if ispc
    xmmDataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
    c4DataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
    omniH5FileStr = 'G:\My Drive\Research\Projects\Data\omni.h5';
elseif isunix
    xmmDataFolder = [initialize_root_path,'LargeFiles/issiteam2017/'];
    c4DataFolder = [initialize_root_path,'LargeFiles/issiteam2017/'];
    omniH5FileStr = [initialize_root_path,'LargeFiles/omni/omni.h5'];
end

xmmOutputFileStr = [dataPathStr,filesep,xmmOutputFileStr];
c4OutputFileStr = [dataPathStr,filesep,c4OutputFileStr];
%% Inputs
% Cropping time
% timeMinXMM = datenum('2000-02-16 21:16:00');
% timeMaxXMM = datenum('2012-08-30 09:29:00');
% timeMinC4 = datenum('2001-01-09 15:23:00');
% timeMaxC4 = datenum('2015-04-30 23:59:00');
% 
timeMinXMM = datenum('2001-01-09 15:23:00');
timeMaxXMM = datenum('2001-01-10 09:29:00');
timeMinC4 = datenum('2001-01-09 15:23:00');
timeMaxC4 = datenum('2001-01-10 09:29:00');


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

% c4.maginput = zeros(length(c4.time),25);
% [spacecraft.c4.maginput,spacecraft.c4.time] = crop_time(c4.maginput,c4.time,timeMinC4,timeMaxC4);
% 
% xmm.maginput = zeros(length(xmm.time),25);
% [spacecraft.xmm.maginput,spacecraft.xmm.time] = crop_time(xmm.maginput,xmm.time,timeMinXMM,timeMaxXMM);

disp('3/7 Done converting coordinates from GSE to GEO...');
%% OMNI Maginput
[omnimaginput,omniTime] = generate_maginput(omniH5FileStr, timeMinXMM, timeMaxC4);
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
out = generate_foot_point(magFieldModel,100,...
        sysaxes,spacecraft.xmm.time,spacecraft.xmm.GEO(:,1),...
        spacecraft.xmm.GEO(:,2),spacecraft.xmm.GEO(:,3),...
        stopAlt,spacecraft.xmm.maginput);
spacecraft.xmm.footType = out.foot;   
spacecraft.xmm.BGEO = out.BGEO;
spacecraft.xmm.eqGEO = out.eqGEO;
spacecraft.xmm.eqBGEO = out.eqBGEO;
spacecraft.xmm.GDZNorth = onera_desp_lib_rotate(out.footNGEO,'geo2gdz',spacecraft.xmm.time);
spacecraft.xmm.GDZSouth = onera_desp_lib_rotate(out.footSGEO,'geo2gdz',spacecraft.xmm.time);
spacecraft.xmm.BGSE = onera_desp_lib_rotate(spacecraft.xmm.BGEO,'geo2gse',spacecraft.xmm.time);
spacecraft.xmm.eqGSE = onera_desp_lib_rotate(spacecraft.xmm.eqGEO,'geo2gse',spacecraft.xmm.time);
spacecraft.xmm.eqBGSE = onera_desp_lib_rotate(spacecraft.xmm.eqBGEO,'geo2gse',spacecraft.xmm.time);
%%
fileID = fopen(xmmOutputFileStr,'w');
fprintf(fileID,'%s\n','DateTime FootType Bx_GSE By_GSE Bz_GSE Eq_x_GSE Eq_y_GSE Eq_z_GSE Beq_x_GSE Beq_y_GSE Beq_z_GSE');
kipsFormat = "%19s %1u %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f %6.2f %6.2f %6.2f\n";

for i=1:1:length(spacecraft.xmm.time)
    fprintf(fileID,kipsFormat,...
        (datestr(spacecraft.xmm.time(i),'yyyy-mm-dd HH:MM:ss'))',...
        spacecraft.xmm.footType(i),...    
        spacecraft.xmm.BGSE(i,:),...
        spacecraft.xmm.eqGSE(i,:),...
        spacecraft.xmm.eqBGSE(i,:));
end    
fclose(fileID);

disp('Done writing XMM magnetic coordinates output data file...');

disp('5/7 Done calculating XMM...');
%% Cluster Calculation

out = generate_foot_point(magFieldModel,100,...
        sysaxes,spacecraft.c4.time,spacecraft.c4.GEO(:,1),...
        spacecraft.c4.GEO(:,2),spacecraft.c4.GEO(:,3),...
        stopAlt,spacecraft.c4.maginput);
   
spacecraft.c4.footType = out.foot;
spacecraft.c4.BGEO = out.BGEO;
spacecraft.c4.eqGEO = out.eqGEO;
spacecraft.c4.eqBGEO = out.eqBGEO;
spacecraft.c4.GDZNorth = onera_desp_lib_rotate(out.footNGEO,'geo2gdz',spacecraft.c4.time);
spacecraft.c4.GDZSouth = onera_desp_lib_rotate(out.footSGEO,'geo2gdz',spacecraft.c4.time);
spacecraft.c4.BGSE = onera_desp_lib_rotate(spacecraft.c4.BGEO,'geo2gse',spacecraft.c4.time);
spacecraft.c4.eqGSE = onera_desp_lib_rotate(spacecraft.c4.eqGEO,'geo2gse',spacecraft.c4.time);
spacecraft.c4.eqBGSE = onera_desp_lib_rotate(spacecraft.c4.eqBGEO,'geo2gse',spacecraft.c4.time);

%%
fileID = fopen(c4OutputFileStr,'w');
fprintf(fileID,'%s \n','DateTime FootType Bx_GSE By_GSE Bz_GSE Eq_x_GSE Eq_y_GSE Eq_z_GSE Beq_x_GSE Beq_y_GSE Beq_z_GSE');

for i=1:1:length(spacecraft.c4.time)
        fprintf(fileID,kipsFormat,...
        (datestr(spacecraft.c4.time(i),'yyyy-mm-dd HH:MM:ss'))',...
        spacecraft.c4.footType(i),...    
        spacecraft.c4.BGSE(i,:),...
        spacecraft.c4.eqGSE(i,:),...
        spacecraft.c4.eqBGSE(i,:));
    
end    
fclose(fileID);
disp('Done writing C4 magnetic coordinates output data file...');

disp('6/7 Done calculating Cluster...');
disp('7/7 Process Complete');

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
