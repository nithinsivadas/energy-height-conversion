function GW=get_tsyganenko_GW_1(yyyy,programDir,omniASCDir)
%% get_tsyganenko_GW_1.m Runs the fortran code MagParametersProgramONE
%                    to calculate G1-G3, and W1-W6: input parameters to the
%                    tsyganenko magnetic field models 2001, 2003.
%                    Uses .dat/.asc files generated from the code.
%--------------------------------------------------------------------------
% Input
%------
% yyyy -     [yyyy] An integer indicating the year for which you would like
%             to calculate G & W.
%
% programDir - A string indicating the directory where the program
%              MagmodelinputONE.exe is. Which is compiled using g77 from
%              MagmodelinputONE.f manually.
%              Default: '~\github\energy-height-conversion\Tools\...
%              External Tools\Tsyganenko_Parameters\MagParameterProgram-rsw\'
% omniASCDir - The director where the ASC omni files will be downloaded by
%              download_omni_cdf(...,'asc');
%              Default: '~\githun\LargeFiles\omni\ASC\'
%--------------------------------------------------------------------------
% Output
%-------
% GW - A struct array containing the following parameters
%      datenum,ByIMF,BzIMF,Velocity_SW,Density_P,Pressure_dynamic,G:[G1 G2 G3]
%      Status8, kp,akp3,dst,Bz1_6:[Bz1,Bz2,Bz3,Bz4,Bz5,Bz6],W:[W1 W2 W3 W4 W5 W6]
%      Status6
%--------------------------------------------------------------------------
% Modified: 30th May 2019; 28th Nov 2020 - change in omni2 format
% Created : 6th Feb 2017
% Author  : Nithin Sivadas
% Ref     :
% Notes   : The structure filed names don't have the standard variable name
% convention.
% Needs g77 installed, and the External Tool:MagParametersProgram from
% http://virbo.org/svn/virbo/qindenton/MagParameterProgram-rsw/
%--------------------------------------------------------------------------
f=filesep;
if nargin < 3
    omniASCDir = [initialize_root_path,'LargeFiles',f,'omni',f,'ASC',f];
end
if nargin <2 || isempty(programDir)
    programDir = [initialize_root_path,...
        'energy-height-conversion',f,'Tools',f,'External Tools',f,'Tsyganenko_Parameters',f,'MagParameterProgram-rsw',f];
end

% Run fortran application only if GW paramters aren't already stored
%     h = waitbar(12/75,'Calculating G and W values');
    formatSpecOmni2 = '%4d %4d %3d %5d %3d %3d %4d %4d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %9.0f %6.1f %6.0f %6.1f %6.1f %6.3f %6.2f %9.0f %6.1f %6.0f %6.1f %6.1f %6.3f %7.2f %7.2f %6.1f %3d %4d %6d %5d %10.2f %9.2f %9.2f %9.2f %9.2f %9.2f %3d %4d %6.1f %6.1f %6d %6d %5.1f %9.6f %7.4f';
    % Check: https://omniweb.gsfc.nasa.gov/html/ow_data.html#3
    ascOmni2FileStr{1} = ['omni2_',num2str(yyyy-1),'.dat'];
    ascOmni2FileStr{2} = ['omni2_',num2str(yyyy),'.dat'];
    inputFileID1 = fopen([omniASCDir,ascOmni2FileStr{1}],'r');
    inputFileID2 = fopen([omniASCDir,ascOmni2FileStr{2}],'r');
    omni2PrevYr = textscan(inputFileID1,formatSpecOmni2);
    omni2CurrYr = textscan(inputFileID2,formatSpecOmni2);
    fclose(inputFileID1);
    fclose(inputFileID2);

    DOY = [(omni2PrevYr{1,2});(omni2CurrYr{1,2})];
    YEAR = [(omni2PrevYr{1,1});(omni2CurrYr{1,1})];
    HR = [(omni2PrevYr{1,3});(omni2CurrYr{1,3})];
    Kp = [(omni2PrevYr{1,39});(omni2CurrYr{1,39})];
    Dst = [(omni2PrevYr{1,41});(omni2CurrYr{1,41})];
%     waitbar(30/75,h,'Calculating G and W values');
    % Columns: YEAR DOY HR MIN  Kp DST
    matrixOutput = [YEAR, DOY, HR, Kp, Dst];
    % Writing kpdst.lst file
    if ~isfolder([programDir,'1min'])
        if ~ispc
            system(['mkdir "',programDir,'1min"']);
        else
            system(['md "',programDir,'1min"']);
        end
    end
    outputFileID = fopen([programDir,'1min',f,'kpdst.lst'],'w');
    formatSpec = '%4u %3u %2u %2u %5d\n';
    fprintf(outputFileID,formatSpec,matrixOutput');
    fclose(outputFileID);
%     waitbar(49/75,h,'Calculating G and W values');

    %% ASC Files
    ascFileStr{1} = ['omni_min',num2str(yyyy-1),'.asc'];
    ascFileStr{2} = ['omni_min',num2str(yyyy),'.asc'];
    if ~isfile([omniASCDir,ascFileStr{1}]) || ~isfile([omniASCDir,ascFileStr{2}])
        warning('Source ASC files not present. Downloading them...');
        download_omni_cdf(['01 Jan ',num2str(yyyy)],[omniASCDir,'..',f,'..',f],'asc');
    end
%     waitbar(57/75,h,'Calculating G and W values');
    outputFileID = fopen([programDir,'1min',f,'omni_min.asc'],'w');
    inputFileID1 = fopen([omniASCDir,ascFileStr{1}],'r');
    inputFileID2 = fopen([omniASCDir,ascFileStr{2}],'r');

    tempStr = fread(inputFileID1,'*char');
    tempStr1 = fread(inputFileID2,'*char');
    fwrite(outputFileID,tempStr);
    fwrite(outputFileID,tempStr1);
    fclose(outputFileID);
    fclose(inputFileID1);
    fclose(inputFileID2);

    %% Running the MagparametersONE
%     waitbar(60/75,h,'Calculating G and W values: Running fortran code');
    currentPath = pwd;
    cd(programDir);
    disp(num2str(yyyy));
    if isunix
        if ~isfile('MagmodelinputONElinux')
        % If exe file does not exists, then compile
        dos('gfortran -fbackslash -o MagmodelinputONElinux MagmodelinputONE_NS.f');
        % Use MagmidelinputONE_NS.f which is modified to work with gfortran
        end
        dos('./MagmodelinputONElinux');
    else
        dos('MagmodelinputONE');
    end
    cd(currentPath);

    %% Converting .d files into matlab variables
%     waitbar(60/75,h,'Storing G and W values as .mat files');
    dataFileID = fopen([programDir,'1min',f,'WGparametersmin.d'],'r');
    formatSpec = '%s';
    nColumns = 1;
    C_text = textscan(dataFileID, formatSpec, nColumns, 'Delimiter', '\n');
    dataFileFormat = '%u %u %u %u %f %f %f %f %f %f %f %f %u %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f %u';
    dataGW = textscan(dataFileID, dataFileFormat);
    fclose(dataFileID);
    YYYY = double(dataGW{1,1});
    DOY = double(dataGW{1,2});
    Hour = double(dataGW{1,3});
    Min = double(dataGW{1,4});
    tempTime=(datenum(datetime(YYYY,1,0))+DOY+Hour/24 +Min/(60*24));
    tempIndex = 1:length(tempTime);
    timeMin = datenum(datetime(yyyy,1,1));
    timeMax = datenum(datetime(yyyy+1,1,1))-1/(24*60);
    timeIndex = crop_time(tempIndex',tempTime,timeMin,timeMax);
    GW.time = tempTime(timeIndex);
    GW.ByIMF = double(dataGW{1,5}(timeIndex));
    GW.BzIMF = double(dataGW{1,6}(timeIndex));
    GW.Velocity_SW = double(dataGW{1,7}(timeIndex));
    GW.Density_P = double(dataGW{1,8}(timeIndex));
    GW.Pressure_dynamic = double(dataGW{1,9}(timeIndex));
    GW.G = [double(dataGW{1,10}(timeIndex)),double(dataGW{1,11}(timeIndex)),double(dataGW{1,12}(timeIndex))];
    GW.Status8 = double(dataGW{1,13}(timeIndex));
    GW.kp = double(dataGW{1,14}(timeIndex));
    GW.akp3 = double(dataGW{1,15}(timeIndex));
    GW.dst = double(dataGW{1,16}(timeIndex));
    GW.Bz1_6 = [double(dataGW{1,17}(timeIndex)),double(dataGW{1,18}(timeIndex)),...
        double(dataGW{1,19}(timeIndex)),double(dataGW{1,20}(timeIndex)),...
        double(dataGW{1,21}(timeIndex)),double(dataGW{1,22}(timeIndex))];
    GW.W = [double(dataGW{1,23}(timeIndex)),double(dataGW{1,24}(timeIndex)),...
        double(dataGW{1,25}(timeIndex)),double(dataGW{1,26}(timeIndex)),...
        double(dataGW{1,27}(timeIndex)),double(dataGW{1,28}(timeIndex))];
    GW.Status6 = double(dataGW{1,29}(timeIndex));
    GW.fieldNames = C_text;

end
