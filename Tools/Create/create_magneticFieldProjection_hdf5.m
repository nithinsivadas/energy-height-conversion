function status = create_magneticFieldProjection_hdf5(...
    magFieldModel,inputH5FileStr,omniH5FileStr,varargin)
%create_magneticFieldProjection_hdf5 Calculates equatorial magnetic field
%parameters of the magnetically conjugate ionospheric point represented by
%a pixel in optical DASC image or the PFISR energy spectra. 

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addRequired(p,'magFieldModel',validScalarPosNum);
addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));
addRequired(p,'omniH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));

addParameter(p,'options',[0,0,0,0,0]);
addParameter(p,'inputAxesNo',0);
addParameter(p,'pixels',32);

parse(p,magFieldModel,inputH5FileStr,omniH5FileStr,varargin{:});
% Settings
options = p.Results.options;
sysaxes = p.Results.inputAxesNo; % GDZ Coordinates as Input 
nPixels0 = p.Results.pixels; % Number of pixels to generate iso-magnetic-field lines in the ionosphere

% Checking if DASC data is available in the energyFlux hdf5 input file
info = h5info(inputH5FileStr); 
nSections = length(info.Groups);
flag = 0;
multiWaitbar('Creating magneticMap HDF5 section',0);
for i = 1:1:nSections
    temp=info.Groups(i).Name;
    flag=flag+strcmp(temp,'/DASC');
end

% If DASC is available
if flag>0
    thisSectionAddress = '/DASC/';
    lat = readh5_variable_at_time(inputH5FileStr,'lat',thisSectionAddress,[]);
    lon = readh5_variable_at_time(inputH5FileStr,'lon',thisSectionAddress,[]);
    alt = readh5_variable_at_time(inputH5FileStr,'alt',thisSectionAddress,[])/1000; % km
    time = unix_to_matlab_time(h5read(inputH5FileStr,'/DASC/time'));
    latNew = modify_matrix_size(lat,nPixels0,nPixels0);
    lonNew = modify_matrix_size(lon,nPixels0,nPixels0);
    altNew = modify_matrix_size(alt,nPixels0,nPixels0);
% If DASC data is not available use PFISR inverted energy flux data
else
    thisSectionAddress = '/magneticFieldAlignedCoordinates/';
    magGeodeticLatLonAlt = readh5_variable_at_time(inputH5FileStr,...
        'magGeodeticLatLonAlt',thisSectionAddress,[]);
    projectionAlt = readh5_variable_at_time(inputH5FileStr,...
        'projectionAlt',thisSectionAddress,[]); 
    indx=find_altitude(magGeodeticLatLonAlt(:,1,3),projectionAlt);
    lat = magGeodeticLatLonAlt(indx,:,1);
    lon = magGeodeticLatLonAlt(indx,:,2);
    alt = magGeodeticLatLonAlt(indx,:,3);
    time = h5read(inputH5FileStr,'/energyFluxFromMaxEnt/time');
    latNew = interp1(1:length(lat),lat,linspace(1,length(lat),nPixels0));
    lonNew = interp1(1:length(lat),lon,linspace(1,length(lat),nPixels0));
    altNew = interp1(1:length(lat),alt,linspace(1,length(lat),nPixels0));
end

% Storing the coordinates of non-NAN pixels
nTime = length(time);

[maginput,timeMaginput] = generate_maginput(omniH5FileStr,time(1),time(end));
if size(timeMaginput)>1
    maginputInterpolated = interp1(timeMaginput,maginput,time,'nearest','extrap');
else
    maginputInterpolated = maginput;
end
maginputInterpolated = filter_irbem_maginput(magFieldModel,maginputInterpolated);

GDZ(2,:)=latNew(~isnan(latNew));
GDZ(3,:)=lonNew(~isnan(latNew));
GDZ(1,:)=altNew(~isnan(latNew));

nPixels = length(GDZ(1,:));
magFieldModelStr=find_irbem_magFieldModelStr(magFieldModel);
data=get_hdf5_dataformat(nTime,nPixels,magFieldModelStr);

%% Create the hdf5 variables
for i = 1:1:length(data(1,:))
    switch length(data{2,i})
        case 4
            dataSize = fliplr(data{2,i});
            chunkSize = [3, 3, 100, 1];
        case 3
            dataSize = fliplr(data{2,i});
            chunkSize = [3, 100, 1];
        case 2
            dataSize = fliplr(data{2,i});
            chunkSize = [1, 100];
    end
    try 
        h5create(inputH5FileStr,data{1,i},dataSize,'ChunkSize', chunkSize,'Deflate', 9); 
        h5writeatt(inputH5FileStr,data{1,i},'Dimensions',data{3,i});
        h5writeatt(inputH5FileStr,data{1,i},'Descriptions',data{4,i});
        h5writeatt(inputH5FileStr,data{1,i},'Units',data{5,i});
    catch ME; 
    end
end

dkTime = 800; % To allow for low memory
nkTime = ceil(nTime/dkTime);
% dkTime = 2;
% nkTime = 2;
start{1} = 1; start{2} = [1 1]; start{3} = [1 1 1]; start{4} = [1 1 1 1]; 
multiWaitbar('Creating magneticMap HDF5 section',0.2);
multiWaitbar('Calculating Part 1',0);
%% Write the variables in
% Write the time-independent variable in
h5write(inputH5FileStr,data{1,2},GDZ);
h5write(inputH5FileStr,data{1,7},posixtime(datetime(time,'ConvertFrom','datenum')));
% Write the time-dependent variable in
for kTime = 1:nkTime
    
    timeEndIndx = min(kTime*dkTime, nTime);
    timeStartIndx = 1 + (kTime-1)*dkTime;
    nTimeMag = length(timeStartIndx:timeEndIndx);
    count{1} = 1; count{2} = [nPixels nTimeMag]; 
    count{3} = [3 nPixels nTimeMag]; count{4} = [3 3 nPixels nTimeMag];
    itimeTemp = 1;
    
    magEqPointGEO = zeros(nTimeMag, nPixels, 3);
    Bgeo = magEqPointGEO;
    B = zeros(nTimeMag, nPixels);
    gradBmag = magEqPointGEO;
    diffB = zeros(nTimeMag, nPixels, 3, 3);
    Lm = zeros(nTimeMag, nPixels);
    Lstar = zeros(nTimeMag, nPixels);
    Blocal = zeros(nTimeMag, nPixels);
    MLT = zeros(nTimeMag, nPixels);
    
    for itime = timeStartIndx:1:timeEndIndx
        
        % Calculating the coordinates of the conjugate magnetic equatorial coordinate
        [~,magEqPointGEOTemp(:,:)]=onera_desp_lib_find_magequator(magFieldModel,options,...
            sysaxes,time(itime),GDZ(1,:),GDZ(2,:),GDZ(3,:),maginputInterpolated(itime,:));

        % Calculating the Bx,By,Bz,Bmag, Gradient B, and Diff B. 
        [BgeoTemp,BTemp,gradBmagTemp,diffBTemp] = onera_desp_lib_get_bderivs(magFieldModel,options,...
        1,time(itime),magEqPointGEOTemp(:,1),magEqPointGEOTemp(:,2),magEqPointGEOTemp(:,3),...
        maginputInterpolated(itime,:)); 
        
        [LmTemp,LstarTemp,BlocalTemp,~,~,MLTTemp] = onera_desp_lib_make_lstar(magFieldModel,...
            options,sysaxes,time(itime),GDZ(1,:),GDZ(2,:),GDZ(3,:),maginputInterpolated(itime,:)); 
        % Storing these measurments in a media file 

        magEqPointGEO(itimeTemp,:,:) = magEqPointGEOTemp;
        Bgeo(itimeTemp,:,:) = BgeoTemp;
        B(itimeTemp,:) = BTemp;
        gradBmag(itimeTemp,:,:) = gradBmagTemp;
        diffB(itimeTemp,:,:,:) = diffBTemp;
        Lm(itimeTemp,:) = LmTemp;
        Lstar(itimeTemp,:) = LstarTemp;
        Blocal(itimeTemp,:) = BlocalTemp;
        MLT(itimeTemp,:) = MLTTemp;
        itimeTemp = itimeTemp + 1;
        
        multiWaitbar('Calculating Part 2','Increment',1/nTimeMag);
    end
    
    % Write the variables in
    h5write(inputH5FileStr,data{1,1},...
        permute(magEqPointGEO,[3 2 1]),start{3},count{3});
    h5write(inputH5FileStr,data{1,3},...
        permute(Bgeo,[3 2 1]), start{3}, count{3});
    h5write(inputH5FileStr,data{1,4},...
        B',start{2},count{2});
    h5write(inputH5FileStr,data{1,5},...
        permute(gradBmag,[3 2 1]), start{3}, count{3});
    h5write(inputH5FileStr,data{1,6},...
        permute(diffB,[4 3 2 1]), start{4}, count{4});
    
    h5write(inputH5FileStr,data{1,8},Lm',start{2},count{2});
    h5write(inputH5FileStr,data{1,9},Lstar',start{2},count{2});
    h5write(inputH5FileStr,data{1,10},Blocal',start{2},count{2});
    h5write(inputH5FileStr,data{1,11},MLT',start{2},count{2});
    
        start{1} = 1; 
        start{2} = [1 1] + [0 nTimeMag]; 
        start{3} = [1 1 1] + [0 0 nTimeMag]; 
        start{4} = [1 1 1 1] + [0 0 0 nTimeMag];  
        multiWaitbar('Calculating Part 1','Increment',1/nkTime);
end
multiWaitbar('Calculating Part 1',1);
multiWaitbar('Calculating Part 2',1);
multiWaitbar('Creating magneticMap HDF5 section',1);
status = 'Success';
end

function data = get_hdf5_dataformat(nTime, nPixels, magFieldModelStr)
    % get_hdf5_dataformat.m Defines the data format for HDF5 group /magneticMap/
    % for the current function. 
    
    data{1,1} = ['/magneticMap/',magFieldModelStr,'/magEqCoordGEO']; % HDF5 address
    data{2,1} = [nTime, nPixels, 3]; % Size
    data{3,1} = 'nTime x nPixels x nCoords(GEO)'; % Dimensions 
    data{4,1} = 'Conjugate equatorial magnetic-field coordinate in GEO(X,Y,Z)'; 
    data{5,1} = 'X,Y,Z in RE'; 

    data{1,2} = ['/magneticMap/',magFieldModelStr,'/ionosphereCoordGDZ']; % HDF5 address
    data{2,2} = [nPixels, 3]; % Size
    data{3,2} = 'nPixels x nCoords(GDZ)'; % Dimensions 
    data{4,2} = 'Input ionosphere coordinates of pixels in GDZ'; 
    data{5,2} = 'lat,lon,alt in [deg,deg,km]'; 

    data{1,3} = ['/magneticMap/',magFieldModelStr,'/BgeoEq']; % HDF5 address
    data{2,3} = [nTime, nPixels, 3]; % Size
    data{3,3} = 'nTime x nPixels x nDirections(GEO)'; % Dimensions 
    data{4,3} = 'Conjugate equatorial magnetic-field vector in GEO(X,Y,Z)'; 
    data{5,3} = 'Bx,By,Bz in nT'; 

    data{1,4} = ['/magneticMap/',magFieldModelStr,'/BmagEq']; % HDF5 address
    data{2,4} = [nTime, nPixels]; % Size
    data{3,4} = 'nTime x nPixels'; % Dimensions 
    data{4,4} = 'Conjugate equatorial magnetic-field magnitude'; 
    data{5,4} = '|B| in nT';
    
    data{1,5} = ['/magneticMap/',magFieldModelStr,'/gradBmagEq']; % HDF5 address
    data{2,5} = [nTime, nPixels, 3]; % Size
    data{3,5} = 'nTime x nPixels x nGradient(GEO)'; % Dimensions 
    data{4,5} = 'Gradient of magnetide of conjugate equatorial magnetic-field in GEO(dX,dY,dZ)'; 
    data{5,5} = 'dB/dx,dB/dy,dB/dz in nT';
    
    data{1,6} = ['/magneticMap/',magFieldModelStr,'/diffBEq']; % HDF5 address
    data{2,6} = [nTime, nPixels, 3, 3]; % Size
    data{3,6} = 'nTime x nPixels x nDirections(i) x nGradient(j)'; % Dimensions 
    data{4,6} = 'Derivatives of conjugate equatorial magnetic-field vector'; 
    data{5,6} = 'dBi/dBj in nT, i&j in GEO';
    
    data{1,7} = ['/magneticMap/',magFieldModelStr,'/time']; % HDF5 address
    data{2,7} = [nTime, 1]; % Size
    data{3,7} = 'nTime x 1'; % Dimensions 
    data{4,7} = 'Time instances'; 
    data{5,7} = 'POSIX';
    
    data{1,8} = ['/magneticMap/',magFieldModelStr,'/Lm']; % HDF5 address
    data{2,8} = [nTime, nPixels]; % Size
    data{3,8} = 'nTime x nPixels'; % Dimensions 
    data{4,8} = 'L McIlwain '; 
    data{5,8} = 'L in [RE]';
    
    data{1,9} = ['/magneticMap/',magFieldModelStr,'/Lstar']; % HDF5 address
    data{2,9} = [nTime, nPixels]; % Size
    data{3,9} = 'nTime x nPixels'; % Dimensions 
    data{4,9} = 'L Roederer  or ?=2?*Bo*/Lstar [nT Re2] '; 
    data{5,9} = 'L* in [RE]';
    
    data{1,10} = ['/magneticMap/',magFieldModelStr,'/BIonosphere']; % HDF5 address
    data{2,10} = [nTime, nPixels]; % Size
    data{3,10} = 'nTime x nPixels'; % Dimensions 
    data{4,10} = 'Magnitude of Bfield at the ionosphere '; 
    data{5,10} = '[nT]';
    
    data{1,11} = ['/magneticMap/',magFieldModelStr,'/MLT']; % HDF5 address
    data{2,11} = [nTime, nPixels]; % Size
    data{3,11} = 'nTime x nPixels'; % Dimensions 
    data{4,11} = 'MLT of the Coordinate'; 
    data{5,11} = '[Hr]';
end




































































