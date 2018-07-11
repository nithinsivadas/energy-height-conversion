%%
clear variables
%%
c4DataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
c4MatFile = [c4DataFolder,'Mat\cluster_all_records.mat'];
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
%%
xmmDataFolder = 'G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';

xmmMatFile = [xmmDataFolder,'Mat\xmm_all_records.mat'];
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
%%


% Cropping time
% timeMin = datenum('2001-01-09 15:23:00');
% timeMax = datenum('2001-01-09 15:27:00');
timeMin = datenum('2000-02-16 21:16:00');
timeMax = datenum('2000-02-16 21:49:00');
[spacecraft.c4.xGSE,spacecraft.c4.time] = crop_time(c4.xGSE,c4.time,timeMin,timeMax);
[spacecraft.xmm.xGSE,spacecraft.xmm.time] = crop_time(xmm.xGSE,xmm.time,timeMin,timeMax);
time = spacecraft.xmm.time;

% Conversion from GSE to GEO
c=define_universal_constants;
% RE = c.RE/1000; %in km
spacecraft.c4.GEO   = onera_desp_lib_rotate(spacecraft.c4.xGSE,'gse2geo',spacecraft.c4.time);
spacecraft.xmm.GEO   = onera_desp_lib_rotate(spacecraft.xmm.xGSE,'gse2geo',spacecraft.xmm.time);

% Calculating maginput/OMNI data
% nMonth = months(datestr(time(1),'mmm dd yyyy'),datestr(time(end),'mmm dd yyyy'))+1;
% dateVec = datevec(time(1));
% dateVec(3) = 1; %setting day to 1st
% omniDataProcessTimes = datenum(dateVec);
% omniDataStorage='G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
% omnimaginput = [];
% omniTime = [];
%     for iMonth = 1:1:nMonth
%         omniDataTemp = process_omni_data(datestr(omniDataProcessTimes,'yyyy mm dd'),omniDataStorage,false);
%         omnimaginput = [omnimaginput; omniDataTemp.minutely.maginput];
%         omniTime = [omniTime; omniDataTemp.minutely.time];
%         omniDataProcessTimes = addtodate(omniDataProcessTimes,iMonth,'month');
%     end
% maginput = interp1(omniTime,omnimaginput,time);
c4.maginput = zeros(length(c4.time),25);
c4.maginput(:,2) =c4.Dst; 
c4.maginput(:,5) = c4.Pdyn;
c4.maginput(:,6) = c4.BYimf; 
c4.maginput(:,7) = c4.BZimf;
[spacecraft.c4.maginput,spacecraft.c4.time] = crop_time(c4.maginput,c4.time,timeMin,timeMax);
% xmm.Pdyn(1:5) = [1.93; 1.16; 1.15; 1.35; nan];
xmm.maginput = zeros(length(xmm.time),25);
xmm.maginput(:,2) =xmm.Dst; 
xmm.maginput(:,5) = xmm.Pdyn;
BimfXMM = onera_desp_lib_rotate([zeros(length(xmm.BYimf),1), xmm.BYimf, xmm.BZimf], 'gse2gsm', xmm.time);
xmm.maginput(:,6) = BimfXMM(:,2); 
xmm.maginput(:,7) = BimfXMM(:,3);
[spacecraft.xmm.maginput,spacecraft.xmm.time] = crop_time(xmm.maginput,xmm.time,timeMin,timeMax);

%% OMNI Maginput
nMonth = months(datestr(timeMin,'mmm dd yyyy'),datestr(timeMax,'mmm dd yyyy'))+1;
dateVec = datevec(timeMin);
dateVec(3) = 1; %setting day to 1st
omniDataProcessTimes = datenum(dateVec);
omniDataStorage='G:\My Drive\Research\Research Trips\2018 April Bern ISSI\Work\Data\';
omnimaginput = [];
omniTime = [];
    for iMonth = 1:1:nMonth
        omniDataTemp = process_omni_data(datestr(omniDataProcessTimes,'yyyy mm dd'),omniDataStorage,true);
        omnimaginput = [omnimaginput; omniDataTemp.minutely.maginput];
        omniTime = [omniTime; omniDataTemp.minutely.time'];
        omniDataProcessTimes = addtodate(omniDataProcessTimes,iMonth,'month');
    end
xmmMaginput = interp1(omniTime,omnimaginput,spacecraft.xmm.time);
c4Maginput = interp1(omniTime,omnimaginput,spacecraft.c4.time);
spacecraft.xmm.maginput = xmmMaginput;
spacecraft.c4.maginput = c4Maginput;
%%

magFieldModel=7;
options = [0,0,0,0,0];
sysaxes = 1; %GEO Input coordinates
hemiFlag = +1;
stopAlt = 110;

spacecraft.c4.BGEO=onera_desp_lib_get_field(magFieldModel,options,...
        sysaxes,spacecraft.c4.time,spacecraft.c4.GEO(:,1),...
        spacecraft.c4.GEO(:,2),spacecraft.c4.GEO(:,3),...
        spacecraft.c4.maginput);
spacecraft.c4.BGSE = onera_desp_lib_rotate(spacecraft.c4.BGEO,'geo2gse',spacecraft.c4.time);

spacecraft.c4.GDZNorth=onera_desp_lib_find_foot_point(magFieldModel,options,...
    sysaxes,spacecraft.c4.time,spacecraft.c4.GEO(:,1),...
    spacecraft.c4.GEO(:,2),spacecraft.c4.GEO(:,3),...
    stopAlt,hemiFlag,spacecraft.c4.maginput);

spacecraft.xmm.BGEO=onera_desp_lib_get_field(magFieldModel,options,...
    sysaxes,spacecraft.xmm.time,spacecraft.xmm.GEO(:,1),...
    spacecraft.xmm.GEO(:,2),spacecraft.xmm.GEO(:,3),...
    spacecraft.xmm.maginput);
spacecraft.xmm.BGSE = onera_desp_lib_rotate(spacecraft.xmm.BGEO,'geo2gse',spacecraft.xmm.time);

spacecraft.xmm.GDZNorth=onera_desp_lib_find_foot_point(magFieldModel,options,...
    sysaxes,spacecraft.xmm.time,spacecraft.xmm.GEO(:,1),...
    spacecraft.xmm.GEO(:,2),spacecraft.xmm.GEO(:,3),...
    stopAlt,hemiFlag,spacecraft.xmm.maginput);

%%
% conjunction.probeStr='xmm';
% conjunction.radius = 500;
% 
% probes = find_magnetic_conjunctions_at_mag_equator(time,spacecraft,conjunction,110,1,magFieldModel);
% probes.time = time;
multiWaitbar('closeall');

% Replace string in file
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
