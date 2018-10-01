function [time,timePrimary] = create_energy_spectra_images(h5FileStr,...
    imageStoreDir,videoFileNameStr,energySlice,eFluxLim, opticalLim, timeMinStr,timeMaxStr,latLim,lonLim,setStoreImage)
    %create_energy_spectra_images.m Create energy spectra images
    %   Creates images of energy spectra overlaid on optical images from
    %   eneryFlux HDF5 file, and then stores the images on a folder specified
    %   by imageStore Dir. 

    if nargin < 11
        setStoreImage = true;
    end
    if nargin < 10 || isempty(lonLim)
        lonLim = [-153 -143];
    end
    if nargin < 9 || isempty(latLim)
        latLim = [63 67];
    end
    if nargin < 8
        timeMaxStr = [];
    end
    if nargin < 7
        timeMinStr = [];
    end
%%
% Generate time array
    % Structure of time variable
    % time->pfisr->h5Address
    %            ->min->indx
    %            ->max->indx
    %      ->dasc->h5Address(nDays)
    %            ->date(nDays) (the date of the DASC group in datenum)
    %            ->min->indx
    %                 ->thisiDay
    %            ->max->indx
    %                 ->thisiDay
    if ~isfolder(imageStoreDir)
        mkdir(imageStoreDir);
    end
    time.pfisr.h5Address = '/energyFluxFromMaxEnt/time';
    time.dasc.h5Address = '/DASC/time';
    thisDascAddress = '/DASC/';
% Read time arrays
    time.pfisr.value(:,1) = h5read(h5FileStr,time.pfisr.h5Address)';
    if isempty(timeMinStr)
        timeMinStr = datestr(min(time.pfisr.value));
    end
    if isempty(timeMaxStr)
        timeMaxStr = datestr(max(time.pfisr.value));
    end
    time.pfisr.min.indx = find_time(time.pfisr.value(:,1),timeMinStr);
    time.pfisr.max.indx = find_time(time.pfisr.value(:,1),timeMaxStr);
    
    time.dasc.value(:,1) = unix_to_matlab_time(h5read(h5FileStr,time.dasc.h5Address)');
  
    time.dasc.min.indx = find_time(time.dasc.value(:,1),timeMinStr);
    time.dasc.max.indx = find_time(time.dasc.value(:,1),timeMaxStr);
%%
% Calculating primary time array  
    if time.pfisr.value(2)-time.pfisr.value(1)<=time.dasc.value(2)-time.dasc.value(1)
        timePrimary = time.pfisr.value(time.pfisr.min.indx:time.pfisr.max.indx);
        whichPrimaryTime = 'pfisr';
    else
        timePrimary = time.dasc.value(time.dasc.min.indx:time.dasc.max.indx);
        whichPrimaryTime = 'dasc';
    end
    nTime = length(timePrimary);
%%
% Reading the important variables
%     dascDayIndxPrimary = floor(timePrimary)-time.dasc.date(1)+1;
    energyBin = readh5_variable_at_time(h5FileStr,'energyBin',...
        '/energyFluxFromMaxEnt/',[])';
    
    magcoords = permute(readh5_variable_at_time(h5FileStr,'magGeodeticLatLonAlt',...
        '/magneticFieldAlignedCoordinates/',[]),[3 2 1]);
    projectionAlt = readh5_variable_at_time(h5FileStr,'projectionAlt',...
        '/magneticFieldAlignedCoordinates/',[]);
    pfisrAltitude = squeeze(magcoords(3,:,:)); 
    % Interpolate beam-wise at projectionAlt, and get the lat, lon values for [nBeamsxnEnergy]
    tempLatitude = squeeze(magcoords(1,:,:)); 
    tempLongitude = squeeze(magcoords(2,:,:));
    
    peakIonizationAlt = calculate_peak_altitude_of_ionization(energySlice*1000,...
        (timePrimary(nTime)+timePrimary(1))/2,tempLatitude(1,:)',tempLongitude(1,:)',...
        pfisrAltitude(1,:)');
    fprintf('PFISR Energy Flux Map projected at ',num2str(peakIonizationAlt),' km');
    nBeams = size(tempLatitude,1);
    for iBeam = 1:1:nBeams
        pfisrLatitude(iBeam,1) = interp1(pfisrAltitude(iBeam,:),tempLatitude(iBeam,:),peakIonizationAlt,'linear','extrap');
        pfisrLongitude(iBeam,1) = interp1(pfisrAltitude(iBeam,:),tempLongitude(iBeam,:),peakIonizationAlt,'linear','extrap');  
    end
    pfisrLatitude = repmat(pfisrLatitude,size(energyBin));
    pfisrLongitude = repmat(pfisrLongitude,size(energyBin));
    multiWaitbar('Store Images',0);
    
    waitBarIncrement = 1./nTime;
    for itime = 1:1:nTime
        thisTime = timePrimary(itime);
        pfisrTimeIndx = find_time(time.pfisr.value,thisTime);
        dascTimeIndx = find_time(time.dasc.value,thisTime);

        pfisrTime = time.pfisr.value(pfisrTimeIndx);
        dascTime = time.dasc.value(dascTimeIndx);
    %     message = readh5_variable_at_time(h5FileStr,'message',thisDascAddress,dascTimeIndx);
        diffEnergyFlux = readh5_variable_at_time(h5FileStr,'energyFlux',...
            '/energyFluxFromMaxEnt/',pfisrTimeIndx)';
        diffEnergyFlux(diffEnergyFlux<0) = 10^3;
        zEnergyBin=repmat(energyBin,size(diffEnergyFlux,1),1);
        
        opticalData = readh5_variable_at_time(h5FileStr,'ASI',thisDascAddress,dascTimeIndx);
        dascLatitude = readh5_variable_at_time(h5FileStr,'lat',thisDascAddress,dascTimeIndx);
        dascLongitude = readh5_variable_at_time(h5FileStr,'lon',thisDascAddress,dascTimeIndx);

        % Plotting
        figureHandle=figure;
        if (abs(dascTime-thisTime)*24*60 > 3) opticalData = NaN(size(opticalData)); end
        if (abs(pfisrTime-thisTime)*24*60 > 3) diffEnergyFlux = NaN(size(diffEnergyFlux)); end
        plot_2D_energy_slice_with_DASC_v2018(figureHandle,thisTime,pfisrTime,...
        diffEnergyFlux, pfisrLatitude, pfisrLongitude,...
            zEnergyBin, energySlice, eFluxLim, opticalData(:), opticalLim, dascTime,...
            dascLatitude(:), dascLongitude(:),...
        imageStoreDir,latLim,lonLim,setStoreImage);
        multiWaitbar('Store Images','Increment',waitBarIncrement);
    end
    
    tempDir = strsplit(imageStoreDir,filesep);
    iImageDir=length(find(~cellfun(@isempty,tempDir)));
    if isunix 
        iImageDir=iImageDir+1; 
    end
    create_video(strjoin(tempDir(1:iImageDir-1),filesep),tempDir{iImageDir},videoFileNameStr);
    multiWaitbar('Store Images','Close');
end

%% Function to plot 2D energy slice overlaid with DASC
% pfisrData.thisTime
% pfisrData.diffEnergyFlux
% pfisrData.latitude
% pfisrData.longitude
% pfisrData.zEnergyBin
% pfisrData.energySlice
% pfisrData.eFluxLim
% dascData.opticalData
% dascData.opticalLim
% dascData.thisTime
% dascData.latitude
% dascData.longitude

function combine_2D_plots_v2018(figureHandle,thisTime,pfisrData,...
        dascData, magFieldData, storeImageDir,latLim,lonLim, setStoreImage)
% Need to write script to provide choice between the several overlays
% Store settings for the plot
settings.energySlice.setTimeLabel = true;

if nargin < 9
    setStoreImage = false;
end
if nargin < 8
    lonLim = [-153 -143];
end
if nargin < 7
    latLim = [63 67];
end
%     eFluxLim = [8 10];
%     eFluxLim = [300 600];

    imageName = strcat('figure_',datestr(thisTime,'HH_MM_SS'));
    if isfile(strcat(storeImageDir,imageName,'.png')) && setStoreImage==true
    warning('Image file already exists, skipping calculation...');
        % Create two axes
    else
    resize_figure(figureHandle,148,210); %A5 Paper Size

    % Plot optical data
    [axesHandleOptical, h1]=plot_DASC_geodetic(dascData.opticalData, dascData.thisTime,...
        dascData.latitude, dascData.longitude, 512, latLim, lonLim);
    % Zoomed coordinates: [64.85 65.05], [-147.95 -147.35]
    % Total : [63 67], [-153 -143]
    
    % Plot magnetic field
    [axesHandleMagnetic, h3] = plot_2D_magnetic_foot_points(magFieldData.magEqPointGEO,...
        magFieldData.ionosphereCoord, magFieldData.BfieldModelStr, magFieldData.thisTime,...
        magFieldData.setMapOn, latLim, lonLim, magFieldData.contourArray,...
        magFieldData.setFieldLabelOn, magFieldData.setTimeLabelOn);
        
    % Plot energy data
    axesHandleEnergy = axes;
    axesm('lambertstd','MapLatLimit',getm(axesHandleOptical,'MapLatLimit'),...
            'MapLonLimit',getm(axesHandleOptical,'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off')
    axis off
    latWidth = latLim(2)-latLim(1);
    lonWidth = lonLim(2)-lonLim(1);
%     pfisrData.diffEnergyFlux(pfisrData.diffEnergyFlux(:)<0)=nan;
    pfisrData.diffEnergyFlux(imag(pfisrData.diffEnergyFlux(:))~=0)=nan;
    [h2]=plot_2D_energy_slice_geodetic_v2018(pfisrData.diffEnergyFlux,...
        pfisrData.latitude, pfisrData.longitude,...
        pfisrData.zEnergyBin, pfisrData.thisTime,...
        pfisrData.energySlice,latWidth,lonWidth,false,settings.energySlice.setTimeLabel);

    colormap(axesHandleOptical,'viridis');
    colormap(axesHandleEnergy,'inferno');
    
    cb1 = colorbar(axesHandleOptical,'eastoutside');
    cb2 = colorbar(axesHandleEnergy,'westoutside');
    ylabel(cb1, '[Rayleigh]');                  
    ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
    caxis(axesHandleOptical,dascData.opticalLim); %clim for optical intensity
    caxis(axesHandleEnergy,pfisrData.eFluxLim);     %clim for energy flux
    % Making Energy spectra translucent
    alpha(axesHandleEnergy,0.5);
    % Linking axes together
    set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
    set(axesHandleMagnetic,'Position',get(axesHandleOptical,'Position'));
    linkaxes([axesHandleOptical,axesHandleEnergy,axesHandleMagnetic]);

% Storing image
    if setStoreImage == true
        export_fig(strcat(storeImageDir,imageName,'.png'),'-r300','-png','-nocrop');
        close(figureHandle);
    end
    end
end



%% Function to plot 2D energy slice overlaid with DASC 
function plot_2D_energy_slice_with_DASC_v2018(figureHandle,thisTime,pfisrTime,...
    diffEnergyFlux, pfisrLatitude, pfisrLongitude,...
        zEnergyBin, energySlice, eFluxLim, opticalData,opticalLim, dascTime,...
        dascLatitude, dascLongitude,...
    storeImageDir,latLim,lonLim, setStoreImage)

% Store settings for the plot
settings.energySlice.setTimeLabel = true;

if nargin < 17
    setStoreImage = false;
end
if nargin < 16
    latLim = [63 67];
end
if nargin < 15
    lonLim = [-153 -143];
end
if isempty(eFluxLim)
    eFluxLim = [8 10];
end
if isempty(opticalLim)
    eFluxLim = [300 600];
end
    imageName = strcat('figure_',datestr(thisTime,'HH_MM_SS'));
    if isfile(strcat(storeImageDir,imageName,'.png')) && setStoreImage==true
    warning('Image file already exists, skipping calculation...');
        % Create two axes
    else
    resize_figure(figureHandle,148,210); %A5 Paper Size

    % Plot optical data
    [axesHandleOptical, h1]=plot_DASC_geodetic(opticalData, dascTime,...
        dascLatitude, dascLongitude, 512, latLim, lonLim);
    % Zoomed coordinates: [64.85 65.05], [-147.95 -147.35]
    % Total : [63 67], [-153 -143]
    
    % Plot energy data
    axesHandleEnergy = axes;
    axesm('lambertstd','MapLatLimit',getm(axesHandleOptical,'MapLatLimit'),...
            'MapLonLimit',getm(axesHandleOptical,'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off')
    axis off
    latWidth = latLim(2)-latLim(1);
    lonWidth = lonLim(2)-lonLim(1);
%     diffEnergyFlux(diffEnergyFlux(:)<0)=nan;
    diffEnergyFlux(imag(diffEnergyFlux(:))~=0)=nan;
    [h2]=plot_2D_energy_slice_geodetic_v2018(diffEnergyFlux, pfisrLatitude, pfisrLongitude,...
        zEnergyBin, pfisrTime,...
        energySlice,latWidth,lonWidth,false,settings.energySlice.setTimeLabel);

    colormap(axesHandleOptical,'viridis');
    colormap(axesHandleEnergy,'inferno');
    
    cb1 = colorbar(axesHandleOptical,'eastoutside');
    cb2 = colorbar(axesHandleEnergy,'westoutside');
    ylabel(cb1, '[Rayleigh]');                  
    ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
    caxis(axesHandleOptical,opticalLim); %clim for optical intensity
    caxis(axesHandleEnergy,eFluxLim);     %clim for energy flux
    % Making Energy spectra translucent
    alpha(axesHandleEnergy,0.5);
    % Linking axes together
    set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
    linkaxes([axesHandleOptical,axesHandleEnergy]);

% Storing image
    if setStoreImage == true
        export_fig(strcat(storeImageDir,imageName,'.png'),'-r300','-png','-nocrop');
        close(figureHandle);
    end
    end
end

%% Plotting energy spectra
function [h2] = plot_2D_energy_slice_geodetic_v2018( diffEnergyFlux,...
    latitude, longitude, zEnergyBin, timeNumPFISR, ...
    thisEnergy, latWidth, lonWidth, setMapOn,...
    setTimeLabelOn, imageSize)
%plot_2D_energy_slice_geodetic.m Plot 2D differential energy flux slices 
%from 4-D PFISR data sets on lat, long map
%--------------------------------------------------------------------------
%Input
%-----
% data          : arranged beam-wise
% -> flux       : differential number flux [nE x nTime]
% -> energyFlux : differential energy flux [nE x nTime]
% -> chi2       : Reduced chi-2 of the maximum entropy regression [1 x nTime]
% -> qInvert    : Production rate from inverted energy flux [nh x nTime]
% -> maxIter    : Maximum number of iterations [1 x nTime]
% -> energyBin  : Energy bin values [nE x 1]
% -> time       : Time array [nTime x 1]
% -> alt        : Altitude array [nh x 1]
% -> A          : Production rate vs. number flux matrix [nh x nE]
% -> qInput     : Production rate derived from electron density [nh x nTime]
% amisrData     : 
% -> site       :latitude, longitude and altitude (PFISR location)
% -> magBeamNo  :the beam number/ID that points along the mag. field line 
% magcoords     : arranged non-beam-wise [nh x 3 x nBeams]
% energyBin     : Energy bin values [nE x 1]
% nBeams        : Total number of beams
% timeNo        : Time number of the energy slice to be plotted
% altitude      : Altitude of projection of the energy slice
% energy        : Energy in keV of the differential energy flux to be plotted
% setMapOn      : True => Map axis on
% setTimeLabelOn: True => The time and energy values are printed on the
%                 plot
%--------------------------------------------------------------------------
%Output
%------
% h2 - plot handle
%
%% 
%----------------------------------------------------------------------------
% Modified: 2nd Feb 2018 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%
if nargin < 11
    imageSize = 512; %Hard coded? based on cal file
end
if nargin < 10
    setTimeLabelOn = true;
end
if nargin < 9
    setMapOn = true;
end
if nargin < 8
    lonWidth = 4;
end
if nargin < 7
    latWidth = 2;
end
%% Generating lat, lon, h, energy coordinates
% diffEnergyFlux(energy,beams) - at a time instant
% lat, lon, energyBin, diffenergyflux - for all data points
%% Generating data slice

F = scatteredInterpolant(latitude(:), longitude(:), zEnergyBin(:), diffEnergyFlux(:),'nearest','none');

latLim = [min(latitude(:)) max(latitude(:))];
lonLim = [min(longitude(:)) max(longitude(:))];

latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);
Vq = F({latq,lonq,thisEnergy*1000});
Vq(Vq<=0)=nan;
if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
        'MapLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',1,'MLineLocation',1);
    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
end

h2=pcolorm(latq,lonq,log10(Vq)); 
set(h2,'EdgeColor','none');

if setTimeLabelOn==true
    hold on;
    textm(latLim(2), lonLim(2)+lonWidth/10, ['PFISR: ', num2str(thisEnergy),' keV'],'color','r');
    textm(latLim(2)-latWidth/20, lonLim(2)+lonWidth/10,...
        [datestr(timeNumPFISR,'HH:MM:SS'),' UT'],'color', 'r');
    hold off;
end

end

