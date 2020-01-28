function [poes, Table, outputFile] = POES_extract_binary_faster(inputFile,outputFile)
%POES_extract_binary.m Extracts .bin POES files using the fortan software provided
% and provides 2 second time resolution of data in the form of a matlab
% structure and table. 
% Details of the Fortran Software and Documentation
% Trying to extract POES data
% 1. Install gcc using MinGW in windows: https://blog.cs.wmich.edu/mingw-setup-guide/
% 2. POES binary files: https://satdat.ngdc.noaa.gov/sem/poes/data/raw/swpc/ (upto 2014)
% 3. POES bindary file extracter for NOAA15+: https://satdat.ngdc.noaa.gov/sem/poes/docs/SoftwareNOAA15+/
% 4. POES binary file extracting documentation: https://satdat.ngdc.noaa.gov/sem/poes/docs/sem2_docs/2003/SEM2v1.2.pdf
% We decided to use gfortran Compiler 
% in readArcSub_v2.f make the following changes to get high resolution
% line 1: program testRead(arcfile,outfile)
% --> line 1: program testRead
% After line 22, include
% --> call getarg(1,arcfile)
% --> call getarg(2,outfile)
% Change line 27: do rec = 1,3
% --> do rec = 1,5400
% 
% After this compile by using the following lines:
% gfortran -c readArcSub_v2.f
% gfortran -c readPOES15Plus.f
% gfortran -o readPOES15Plus.exe readPOES15Plus.o readArcSub_v2.o
% ./readPOES15Plus.exe inputFile.bin outputFile.txt

tempStr = strsplit(inputFile,filesep);
tempStr1 = strsplit(tempStr{end},'.');
outputFileName = [tempStr1{1},'.txt'];
if nargin<2
    outputFile = [initialize_root_path,'LargeFiles',filesep,'POES',filesep,outputFileName];
end

POESToolPath = [initialize_root_path,'energy-height-conversion',filesep,...
    'Tools',filesep,'External Tools',filesep,'NOAA15PlusSoftware',filesep];
exe = ['"',POESToolPath,'readPOES15Plus.exe"'];

cmd = [exe,' ',['"',inputFile,'"'],' ',['"',outputFile,'"']];

% Running the Fortran code
[status,dat] = system(cmd);
if status ~= 0, error(dat), end

% Reading output of the Fortran Code
Table = readtable(outputFile);

% Interpolating time, lat and lon arrays for 2 second time bins
n=size(Table,1);
di = 1./n;
% multiWaitbar('POES Interpolating time bins...',0);
indx = 1:1:height(Table);
sec = mod(indx-1,4)';

time=(datetime(Table{:,'year'},1,Table{:,'doy'},...
    Table{:,'hr'},Table{:,'min'},Table{:,'sec'}+2*sec));

% time=(datenum(strcat(num2str(Table{sec==0,'year'}),' 1 ',num2str(Table{sec==0,'doy'}),...
%      ' ',num2str(Table{sec==0,'hr'}),':',num2str(Table{sec==0,'min'}),':',num2str(Table{sec==0,'sec'}))));
lat = Table{:,'lat'};
lon = Table{:,'lon'};
alt = Table{:,'alt_km_'};

lat(sec~=0) = nan;
lon(sec~=0) = nan;
alt(sec~=0) = nan;

% multiWaitbar('POES Interpolating time bins...','Increment',di);


%%
poes.time = time;
poes.lat = interp_nans(lat);
poes.lon = interp_nans(lon);
poes.alt = interp_nans(alt);
poes.mep0P1 = Table.mep0P1;
poes.mep0P2 = Table.mep0P2;
poes.mep0P3 = Table.mep0P3;
poes.mep0P4 = Table.mep0P4;
poes.mep0P5 = Table.mep0P5;
poes.mep0P6 = Table.mep0P6;
poes.mep0E1 = Table.mep0E1;
poes.mep0E2 = Table.mep0E2;
poes.mep0E3 = Table.mep0E3;
poes.mep90P1 = Table.mep90P1;
poes.mep90P2 = Table.mep90P2;
poes.mep90P3 = Table.mep90P3;
poes.mep90P4 = Table.mep90P4;
poes.mep90P5 = Table.mep90P5;
poes.mep90P6 = Table.mep90P6;
poes.mep90E1 = Table.mep90E1;
poes.mep90E2 = Table.mep90E2;
poes.mep90E3 = Table.mep90E3;
poes.ted01 = Table.ted01;
poes.ted02 = Table.ted02;
poes.ted03 = Table.ted03;
poes.ted04 = Table.ted04;
poes.ted05 = Table.ted05;
poes.ted06 = Table.ted06;
poes.ted07 = Table.ted07;
poes.ted08 = Table.ted08;
poes.tedfx1 = Table.tedfx1;
poes.tedfx2 = Table.tedfx2;
poes.tedfx5 = Table.tedfx5;


end

