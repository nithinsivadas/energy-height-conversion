function [beamNo, latitude, longitude, altitudeGrid] =...
    get_pfisr_beam_coordinates(fileNameStr, msrParameter, altitudeGrid)

%% get_pfisr_beam_coordinates.m 
%
% Calculates the coordinates of measurement parameters varying with
% altitude along different beams of PFISR (Poker Flat Incoherent Radar)
% system. 
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

%----------------------------------------------------------------------------
% Output 
%--------
% beamNo      : The beam ID of all beams in the measurment
% latitude    : Geographic latitude value of each measurement along beams
%               at the altitude points specified
% longitude   : Geographic longitude value of each measurement along beams
%               at the altitude points specified
% altitudeGrid: The default altitude (of beam 1) 
%               or user specified altitude points

%%
%----------------------------------------------------------------------------
% Modified: 14th Oct 2016 
% Created : 14th Oct 2016
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
    coordinateNo=1; % Using slant range to evaluate the number of beams
    nBeams = get_nbeams(pfisrGD, coordinateNo);
    
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
            
            r1=pfisrGD.dataloc((i-1)*nBeams+iBeam,1);
            range(i,iBeam)=r1;
            azimuth(i,iBeam)=pfisrGD.dataloc((i-1)*nBeams+iBeam,2);
            elevation(i,iBeam)=pfisrGD.dataloc((i-1)*nBeams+iBeam,3);
            
%             MsrPar(i,iBeam,1:1:nTime)=...
%             10.^data((i-1)*nBeams+iBeam,1:1:nTime);

            Altitude(i,iBeam)=...
            r1*sin(pfisrGD.dataloc((i-1)*nBeams+iBeam,end)*pi/180);
            
        end;
    
    end;
    
    %% Coordinate transformation to lat, long, altitude
    pfisrloc=[65.12992 -147.47104 213]; % [deg, deg, m]
    for iBeam=1:1:nBeams
        [lat(:,iBeam), lon(:,iBeam), alt(:,iBeam)] = aer2geodetic( ...
    azimuth(:,iBeam), elevation(:,iBeam), range(:,iBeam)*1000, pfisrloc(1), pfisrloc(2), pfisrloc(3),wgs84Ellipsoid);
    end;

% Assigning the default altitude grid (vertical projection of the measurement along magnetic field aligned beam)
     if (setDefaultAltitudeGrid==true) 
         altitudeGrid=Altitude(:,1); %Along the magnetic field line
     end;

% Interpolating the lat long coordinates along altitude points specified by altitudeGrid
     latitude=zeros(size(altitudeGrid,1),nBeams);
     longitude=zeros(size(altitudeGrid,1),nBeams);
     for iBeam=1:1:nBeams
         latitude(:,iBeam,:)=interp1...
         (Altitude(:,iBeam),lat(:,iBeam,:),altitudeGrid,'linear','extrap');
         longitude(:,iBeam,:)=interp1...
         (Altitude(:,iBeam),lon(:,iBeam,:),altitudeGrid,'linear','extrap');
     end;
     
     beamNo=1:1:nBeams;

    [isThereNAN, totalNAN] = check_nan(latitude);
end   



