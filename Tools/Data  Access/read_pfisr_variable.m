function [msrParValue, altitudeGrid, time] =...
    read_pfisr_variable(fileNameStr, msrParameter, mode, altitudeGrid)

%% read_pfisr_variable.m 
% Calculates the value of the parmater specified in 'msrParameter'
% along altitude from measurements along the different beams of PFISR (Poker
% Flat Incoherent Radar) system. 
%
%------------------------------------------------------------------------------
% Input
%-------
% fileNameStr : [String] The path to HDF5 file which contains the PFISR data
%
% msrParameter: [String] - The name of the measurement parameter
%               (specified in the HDF5 file) which you would like to
%               extract from the HDF5 file (e.g. 'dpopl','popl')
%               Default: Uncorrected electron density: 'popl'
%               IMPORTANT- The parameter has to be expressed in log scale
% altitudeGrid: [N x 1] [km] - A column vector of the altitude points for 
%               which you need the 
%               value of the measurement parameter 
%               Default value of the altitude grid is the altitude projection 
%               of the range points of the beam pointing along 
%               the magnetic field line
% mode        : Determines the kind of data output [Integer values]
%         '-1': Retrieves data along the magnetic field aligned beam
%          '0': [Default mode] Retrieves data averaged across the total number 
%               of beams 
%    '1' or >1: Assumes the number as the beam number along which you would 
%               like the data to be retrieved
%----------------------------------------------------------------------------
% Output 
%--------
% msrParValue : The value of the measurement parameter per altitude; 
%               Units as given by the HDF5 file
% altitudeGrid: The default altitude (of beam 1) 
%               or user specified altitude points
% time        : Time array
%%
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

% Issues: 
% 1. It is assumed that the total number of measurement along every beam is the same, and they share the same time array
% 2. Need to initialize GeoData; Use this function initialize_geodata(filePath)
    
    switch nargin
        case 4
            setDefaultAltitudeGrid=false;  
        case 3
            setDefaultAltitudeGrid=true;  
        case 2
            setDefaultAltitudeGrid=true;
            mode=0;
        case 1
            setDefaultAltitudeGrid=true;
            mode=0;
            msrParameter='popl';
        otherwise
            if nargin>4
                error('TooManyInputs: Number of input variables should be at most 4')
            else
                error('TooFewInputs: Number of input variables should at least be 2');
            end
    end


    % Invoking the GeoData class to extract range and measurement parameter from the HDF5 file
    pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range',msrParameter});
    
    % Estimating the time array
    time=unixtime2matlab(pfisrGD.times(:,1));
    
    % Identifying the number of beams
    temp = pfisrGD.dataloc(1,1);
    j=1;
    while (pfisrGD.dataloc(j+1,1)==temp)
        j=j+1;
    end;    
    nBeams=j; % Total number of beams

    % Estimating the length of the range and time array
    nRange = size(pfisrGD.dataloc,1);
    nTime = length(time);

    % Storing the measured parameter values into a variable
    eval(sprintf('data=pfisrGD.data.%s;',(msrParameter)));
    
    % Initializing certain variables
    Altitude=zeros(nRange/nBeams-1,1);
    MsrPar=zeros(nRange/nBeams-1,nBeams,nTime);
        
    initialBeam=1; %Initial beam number

    % Calculating the altitude for each beam
    for iBeam=1:1:nBeams
        
        for i=1:1:(nRange/nBeams - 1)
            
            r1=pfisrGD.dataloc((i-1)*nBeams+initialBeam,1);
            
            MsrPar(i,iBeam,1:1:nTime)=...
            10.^data((i-1)*nBeams+iBeam,1:1:nTime);

            Altitude(i,iBeam)=...
            r1*sin(pfisrGD.dataloc((i-1)*nBeams+iBeam,end)*pi/180);

        end;
    
    end;

    % Assigning the default altitude grid (vertical projection of the measurement along magnetic field aligned beam)
    if (setDefaultAltitudeGrid==true) 
        altitudeGrid=Altitude(:,1); %Along the magnetic field line
    end;
    
    % Interpolating the measurement parameter along altitude points specified by altitudeGrid
    msrParAlongGrid=zeros(size(altitudeGrid,1),nBeams,nTime);
    for iBeam=1:1:nBeams
        msrParAlongGrid(:,iBeam,:)=interp1...
        (Altitude(:,iBeam),MsrPar(:,iBeam,:),altitudeGrid,'linear','extrap');
    end;
    
    % Estimating the Measurement Parameter Value as specified by the mode

    if mode==-1 % Along magnetic field line
        beamNo=1;
        msrParValue=squeeze(msrParAlongGrid(:,beamNo,:));

    elseif mode==0 % Average of all beams, projected along altitude grid
    % Averaging measurement parameter per altitude slice
        
        msrParSumOvrBeams=...
        zeros(size(msrParAlongGrid,1),size(msrParAlongGrid,3));
        msrParMeanOvrBeams=...
        zeros(size(msrParAlongGrid,1),size(msrParAlongGrid,3));
    
        for i=1:1:size(msrParAlongGrid,1)
            msrParSumOvrBeams(i,:)=squeeze(nansum(msrParAlongGrid(i,:,:),2))';
            msrParMeanOvrBeams(i,:)=msrParSumOvrBeams(i,:)./...
            (nBeams*ones(1,nTime)...
                -squeeze(sum(isnan(msrParAlongGrid(i,:,:)),2))');
        end;
        msrParValue=msrParMeanOvrBeams;

    elseif mode>=1 % Along the particular beam number specified by  mode
        beamNo=mode;
        msrParValue=squeeze(msrParAlongGrid(:,beamNo,:));

    end;
    [isThereNAN, totalNAN] = check_nan(msrParValue);
end   



