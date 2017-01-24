function [msrParValue, xEastAlongGrid, yNorthAlongGrid, altitudeGrid, time] =...
    extract_pfisr_beam_data(pfisrGD, msrParameter, beamNo, nBeams, altitudeGrid)

%% Extract PFISR Beam Data, given input coordinates are in xEast, yNorth, zUP

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
        case 5
           setDefaultAltitudeGrid=false; 
        case 4
           setDefaultAltitudeGrid=true; 
        case 3
           setDefaultAltitudeGrid=true; 
           nBeams = 26;
        case 2
            beamNo=1;
            setDefaultAltitudeGrid=true; 
            nBeams = 26;
        case 1
            setDefaultAltitudeGrid=true;
            nBeams = 26;
            beamNo=1;
            msrParameter='popl';
        otherwise
            if nargin>5
                error('TooManyInputs: Number of input variables should be at most 3')
            else
                error('TooFewInputs: Number of input variables should at least be 0');
            end
    end

    % Estimating the time array
    time=unixtime2matlab(pfisrGD.times(:,1));
    
    % Estimating the length of the range and time array
    nRange = size(pfisrGD.dataloc,1);
    nTime = length(time);

    % Storing the measured parameter values into a variable
    eval(sprintf('data=pfisrGD.data.%s;',(msrParameter)));
    
    % Initializing certain variables
    Altitude=zeros(nRange/nBeams-1,1);
    MsrPar=zeros(nRange/nBeams-1,nBeams,nTime);
    xEast=zeros(nRange/nBeams-1,nBeams);
    yNorth=zeros(nRange/nBeams-1,nBeams);
    xEastAlongGrid=zeros(nRange/nBeams-1,nBeams);
    yNorthAlongGrid=zeros(nRange/nBeams-1,nBeams);
        
    initialBeam=1; %Initial beam number

    % Calculating the altitude for each beam
    for iBeam=1:1:nBeams
        
        for i=1:1:(nRange/nBeams - 1)
                           
            MsrPar(i,iBeam,1:1:nTime)=data((i-1)*nBeams+iBeam,1:1:nTime);

            Altitude(i,iBeam)=pfisrGD.dataloc((i-1)*nBeams+iBeam,end);
            xEast(i,iBeam)=pfisrGD.dataloc((i-1)*nBeams+iBeam,1);
            yNorth(i,iBeam)=pfisrGD.dataloc((i-1)*nBeams+iBeam,2);
        end;
    
    end;

    % Assigning the default altitude grid (vertical projection of the measurement along magnetic field aligned beam)
    if (setDefaultAltitudeGrid==true) 
        altitudeGrid=Altitude(:,beamNo); %Along the magnetic field line
    end;
        % Interpolating the measurement parameter along altitude points specified by altitudeGrid
    msrParAlongGrid=zeros(size(altitudeGrid,1),nBeams,nTime);
    for iBeam=1:1:nBeams
        msrParAlongGrid(:,iBeam,:)=interp1...
        (Altitude(:,iBeam),MsrPar(:,iBeam,:),altitudeGrid,'linear','extrap');
        xEastAlongGrid(:,iBeam)=interp1...
        (Altitude(:,iBeam),xEast(:,iBeam),altitudeGrid,'linear','extrap');
        yNorthAlongGrid(:,iBeam)=interp1...
        (Altitude(:,iBeam),yNorth(:,iBeam),altitudeGrid,'linear','extrap');
    end;
        
    msrParValue=squeeze(msrParAlongGrid(:,beamNo,:));

    [isThereNAN, totalNAN] = check_nan(msrParValue);
end   




