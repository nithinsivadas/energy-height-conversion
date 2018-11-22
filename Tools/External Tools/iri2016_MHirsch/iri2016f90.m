function data = iri2016f90(time,altKm,glat,glon,setSamplePlot)
%IRI2016F90 Matlab-Python wrapper of IRI2016 by MHirsch
%   time,glat,glon are one value each, 
%   while altKmRange is [altMin, altMax, StepSize] (3x1 array)
% Issues: Check the consistency with other implimentations. 

% Verifying if variables are of the right type
narginchk(4,5)
validateattributes(altKm, {'numeric'}, {'positive', 'vector'})
validateattributes(glat, {'numeric'}, {'scalar'})
validateattributes(glon, {'numeric'}, {'scalar'})

% Check if a sample plot needs to be made
if nargin<5
    setSamplePlot = false;
end

validateattributes(setSamplePlot, {'logical'}, {'scalar'})

% Convert time into the right format
switch class(time)
    case {'datetime', 'double'}, t = sprintf('%d %d %d %d %d %d',datevec(time));
end


% Generate a uniform altitude grid
Dalt = (2^(nextpow2(min(diff(altKm)))-1)); 
if Dalt<0.0625 Dalt=0.0625; end % Minimum step size is 0.0625
altKmUniform = floor(altKm(1)):Dalt:ceil(altKm(end));
maxRows = 1000; % Maximum rows fortran code can process
Nalt = length(altKmUniform);
Nexe = ceil(Nalt/maxRows);

% Initializing paths to run the fortran IRI2016
rootPath =[initialize_root_path,'IRI2016',filesep','matlab'];
exeStr = ['..',filesep,'bin',filesep,'iri2016_driver'];
currentPath = pwd;
cd(rootPath);

% Chopping up the uniform grid to be within the maximum number of rows
% for each fortran code run
arrUniform = [];    
for i = 1:Nexe
    % Converting altitude array into [minAlt, maxAlt, Step] sizes
    % which the fortran code take as inputs
    altIndx(1) = 1 + (i-1)*maxRows;
    altIndx(2) = min(i*maxRows,Nalt);
    altKmRange(1) = altKmUniform(altIndx(1)); 
    altKmRange(2) = altKmUniform(altIndx(2));
    altKmRange(3) = Dalt;
    
    arrTemp = execute_iri2016_driver(t,altKmRange,glat,glon,exeStr);
    arrUniform = [arrUniform, arrTemp];
end
cd(currentPath); % Go back to the original directory

arr = interp1(altKmUniform',arrUniform',altKm')';

% Store them on a data structure
data.coordinates.altKm = arr(1,:); % ALTITUDE VECTOR [KM]

data.Ne     = arr(2,:); % ELECTRON DENSITY [M-3]
data.Tn     = arr(3,:); % NEUTRAL TEMPERATURE [K]
data.Ti     = arr(4,:); % ION TEMPERATURE [K]
data.Te     = arr(5,:); % ELECTRON TEMPERATURE [K]

data.nOI    = arr(6,:); % O+ ION DENISTY [M-3]
data.nHI    = arr(7,:); % H+ ION DENISTY [M-3]
data.nHeI   = arr(8,:); % He+ ION DENISTY [M-3]
data.nO2I   = arr(9,:); % O2+ ION DENISTY [M-3]
data.nNOI   = arr(10,:);% NO+ ION DENISTY [M-3]
data.nCI    = arr(11,:);% CLUSTER IONS DEN + ION DENISTY [M-3]
data.nNI    = arr(12,:);% N+ ION DENISTY [M-3]

data.totalIonDensity = (data.nOI + data.nCI +...
    data.nHeI + data.nHI + data.nNI + data.nNOI +...
    data.nO2I); % [M-3] Total Ion Number Density

data.time = time;
data.coordinates.lat = glat;
data.coordinates.lon = glon;

data.description = 'SI Units';

% Plot if requested
if setSamplePlot
    plot_iono(data);
end 

end

function arr = execute_iri2016_driver(t,altKmRange,glat,glon,exeStr)
    
    
    % Running Fortran Code
    [status,dat] = system([exeStr, ' ', t,...
                           ' ',num2str(glat),' ',num2str(glon),...
                           ' ',num2str(altKmRange(1)),...
                           ' ',num2str(altKmRange(2)),...
                           ' ',num2str(altKmRange(3))]);
    if status ~= 0, error(dat), end
    
    % Parse the output from the fortran file
    arr = sscanf(dat, '%f', [12,inf]);


end

function plot_iono(data)

figure;
semilogx(data.Ne(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nOI(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nO2I(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nNOI(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nHI(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nHeI(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nCI(:), data.coordinates.altKm(:))
hold on;
semilogx(data.nNI(:), data.coordinates.altKm(:))
set(gca,'xscale','log')
legend('Ne','nO+','nO2+','nNO+','nH+','nHe+','nClusterIons','nN+','Location','eastoutside');
xlim([10^8 10^12]);

title({[datestr(data.time),' deg.  (',num2str(data.coordinates.lat),', ', num2str(data.coordinates.lon),')']})
xlabel('Density [m^-3]')
ylabel('altitude [km]')
grid('on')


end
