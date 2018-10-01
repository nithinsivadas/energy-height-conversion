function status = create_magneticFieldProjection_hdf5(...
    magFieldModel,inputH5FileStr,omniH5FileStr)
%create_magneticFieldProjection_hdf5 Calculates equatorial magnetic field
%parameters of the magnetically conjugate ionospheric point represented by
%a pixel in optical DASC image or the PFISR energy spectra. 

% Settings
options = [0,0,0,0,0];
sysaxes = 0; % GDZ Coordinates as Input 
nPixels = 64; % Number of pixels to generate iso-magnetic-field lines in the ionosphere

% Checking if DASC data is available in the energyFlux hdf5 input file
info = h5info(inputH5FileStr); 
nSections = length(info.Groups);
flag = 0;
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
    latNew = modify_matrix_size(lat,nPixels,nPixels);
    lonNew = modify_matrix_size(lon,nPixels,nPixels);
    altNew = modify_matrix_size(alt,nPixels,nPixels);
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
    latNew = interp1(1:length(lat),lat,linspace(1,length(lat),nPixels));
    lonNew = interp1(1:length(lat),lon,linspace(1,length(lat),nPixels));
    altNew = interp1(1:length(lat),alt,linspace(1,length(lat),nPixels));
end

% Storing the coordinates of non-NAN pixels
nTime = length(time);

[maginput,timeMaginput] = generate_maginput(omniH5FileStr,time(1),time(end));
maginputInterpolated = interp1(timeMaginput,maginput,time,'nearest','extrap');
GDZ(2,:)=latNew(~isnan(latNew));
GDZ(3,:)=lonNew(~isnan(latNew));
GDZ(1,:)=altNew(~isnan(latNew));

nPixels = length(GDZ(1,:));
data=get_hdf5_dataformat(nTime,nPixels);

for itime = 1:1:1
    
    % Calculating the coordinates of the conjugate magnetic equatorial coordinate
    [~,magEqPointGEOTemp(:,:)]=onera_desp_lib_find_magequator(magFieldModel,options,...
        sysaxes,time(itime),GDZ(1,:),GDZ(2,:),GDZ(3,:),maginputInterpolated(itime,:));
    
    % Calculating the Bx,By,Bz,Bmag, Gradient B, and Diff B. 
    [BgeoTemp,BTemp,gradBmagTemp,diffBTemp] = onera_desp_lib_get_bderivs(magFieldModel,options,...
    1,time(itime),magEqPointGEOTemp(:,1),magEqPointGEOTemp(:,2),magEqPointGEOTemp(:,3),...
    maginputInterpolated(itime,:)); 
    
    % Storing these measurments in a media file 
    magEqPointGEO(itime,:,:) = magEqPointGEOTemp;
    Bgeo(itime,:,:) = BgeoTemp;
    B(itime,:,:) = BTemp;
    gradBmag(itime,:,:) = gradBmagTemp;
    diffB(itime,:,:,:) = diffBTemp;
    
% Bgeo: BxGEO,ByGEO and BzGEO of the magnetic field (nT), (array(3,NTIME_MAX) of double) 
% Bmag: magnitude of magnetic field (array(NTIME_MAX) of double) (nT) 
% gradBmag: gradient (GEO) of magnitude of magnetic field (array(3,NTIME_MAX) of double) (nT) 
% diffB: derivatives of magnetic field vector (array(3,3,NTIME_MAX) of double) (nT)
% diffB(i,j,t) = dB_i/dx_j at t'th location. GEO coordinates. 
end

% Create the hdf5 variables
magFieldModelStr=find_irbem_magFieldModelStr(magFieldModel);
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/magEqCoordGEO'],...
        size(permute(magEqPointGEO,[3 2 1]))); catch ME; end
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/ionosphereCoordGDZ'],...
        size(GDZ)); catch ME; end
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/BgeoEq'],...
        size(permute(Bgeo,[3 2 1]))); catch ME; end
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/BmagEq'],...
        size(B')); catch ME; end
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/gradBmagEq'],...
        size(permute(gradBmag,[3 2 1]))); catch ME; end
try h5create(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/diffBEq'],...
        size(permute(diffB,[4 3 2 1]))); catch ME; end
% Write the variables in 
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/magEqCoordGEO'],...
    permute(magEqPointGEO,[3 2 1]));
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/ionosphereCoordGDZ'],...
    GDZ);
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/BgeoEq'],...
    permute(Bgeo,[3 2 1]));
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/BmagEq'],...
    B');
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/gradBmagEq'],...
    permute(gradBmag,[3 2 1]));
h5write(inputH5FileStr,['/magneticMap/',magFieldModelStr,'/diffBEq'],...
    permute(diffB,[4 3 2 1]));

status = 'Success';
end

function data = get_hdf5_dataformat(nTime, nPixels, magFieldModelStr)
    
    data{1,1} = ['/magneticMap/',magFieldModelStr,'/magEqCoordGEO']; % HDF5 address
    data{2,1} = [nTime, nPixels, 3]; % Size
    data{3,1} = 'nTime x nPixels x nCoords(GEO)'; % Dimensions 
    data{4,1} = 'Conjugate equatorial magnetic-field coordinate in GEO(X,Y,Z)'; 
    data{5,1} = 'X,Y,Z in RE'; 

    data{1,2} = ['/magneticMap/',magFieldModelStr,'/ionosphereCoordGDZ']; % HDF5 address
    data{2,2} = [3, nPixels]; % Size
    data{3,2} = 'nCoords(GDZ) x nPixels'; % Dimensions 
    data{4,2} = 'Input ionosphere coordinates of pixels in GDZ'; 
    data{5,2} = 'lat,lon,alt in [deg,deg,km]'; 

    data{1,3} = ['/magneticMap/',magFieldModelStr,'/BgeoEq']; % HDF5 address
    data{2,3} = [nTime, nPixels, 3]; % Size
    data{3,3} = 'nTime x nPixels x nDirections(GEO)'; % Dimensions 
    data{4,3} = 'Conjugate equatorial magnetic-field vector in GEO(X,Y,Z)'; 
    data{5,3} = 'Bx,By,Bz in nT'; 

    data{1,3} = ['/magneticMap/',magFieldModelStr,'/BmagEq']; % HDF5 address
    data{2,3} = [nTime, nPixels]; % Size
    data{3,3} = 'nTime x nPixels'; % Dimensions 
    data{4,3} = 'Conjugate equatorial magnetic-field magnitude'; 
    data{5,3} = '|B| in nT';
    
    data{1,4} = ['/magneticMap/',magFieldModelStr,'/gradBmagEq']; % HDF5 address
    data{2,4} = [nTime, nPixels, 3]; % Size
    data{3,4} = 'nTime x nPixels x nGradient(GEO)'; % Dimensions 
    data{4,4} = 'Gradient of magnetide of conjugate equatorial magnetic-field in GEO(dX,dY,dZ)'; 
    data{5,4} = 'dB/dx,dB/dy,dB/dz in nT';
    
    data{1,5} = ['/magneticMap/',magFieldModelStr,'/diffBEq']; % HDF5 address
    data{2,5} = [nTime, nPixels, 3, 3]; % Size
    data{3,5} = 'nTime x nPixels x nDirections(i) x nGradient(j)'; % Dimensions 
    data{4,5} = 'Derivatives of conjugate equatorial magnetic-field vector'; 
    data{5,5} = 'dBi/dBj in nT, i&j in GEO';
    
    
end







































































