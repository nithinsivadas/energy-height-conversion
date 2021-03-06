function plot_2D_energy_slice_with_DASC(figureHandle,timeStr,...
    dataInvEnergySpectra,amisrData,magcoords,DASCcalFile,DASCdataFolder,...
    latLim,lonLim, storeImageDir)

% Store settings for the plot
settings.energySlice.setTimeLabel = true;

rootPathStr=initialize_root_path;
if nargin < 10
    setStoreImage = false;
else 
    setStoreImage = true;
end
if nargin < 9
    latLim = [63 67];
end
if nargin < 8
    lonLim = [-153 -143];
end
if nargin < 7
    DASCdataFolder = [rootPathStr,'\LargeFiles\DASC\26Mar2008'];
end
if nargin < 6
    DASCcalFile.az = [rootPathStr,'\LargeFiles\DASC\26Mar2008\PKR_Cal_before_2011\PKR_20111006_AZ_10deg.FITS'];
    DASCcalFile.el = [rootPathStr,'\LargeFiles\DASC\26Mar2008\PKR_Cal_before_2011\PKR_20111006_EL_10deg.FITS'];
end

dirName = DASCdataFolder;
dataInv = dataInvEnergySpectra;
azOldRes = fitsread(DASCcalFile.az);
elOldRes = fitsread(DASCcalFile.el);

% Storing file names in the DASC data folder onto fileStr
files = dir(dirName);
fileIndex = find(~[files.isdir]);
fileTempCells  = struct2cell(files);
fileStr = fileTempCells(1,fileIndex);

% Generating Time Stamp
aldtnum = fitsfiletimestamp(fileStr);
timeASI = (aldtnum-datenum('jan-01-1970'))*(24*3600);
timeASI = unix_to_matlab_time(timeASI);

timeIndex = find_time(timeASI,timeStr);

ASIDataStr = strcat(dirName,'\',(fileStr(timeIndex)));

    try
        [ASI.dataNew, ASI.lat, ASI.lon, ASI.az_new, ASI.el_new,...
            ASI.sensorloc, ASI.timeDASC] = DASC_aer_to_geodetic...
            (char(ASIDataStr), azOldRes, elOldRes,...
            512, 30, 110);
        catch ME
          warning(['Could not load file: ',fileStr(timeIndex)]);
    end


 %% Create two axes
    
    resize_figure(figureHandle,148,210); %A5 Paper Size

    % Plot optical data
    [axesHandleOptical, h1]=plot_DASC_geodetic(ASI.dataNew, ASI.timeDASC,...
        ASI.lat, ASI.lon, 512, latLim, lonLim);
    % Zoomed coordinates: [64.85 65.05], [-147.95 -147.35]
    % Total : [63 67], [-153 -143]
    
    % Identify the closest time in the PFISR data to the current optical
    % time
    timePFISRNo = find_time(dataInv(1).time,datestr(timeASI(timeIndex)));
    % Plot energy data
    axesHandleEnergy = axes;
    axesm('lambertstd','MapLatLimit',getm(axesHandleOptical,'MapLatLimit'),...
            'MapLonLimit',getm(axesHandleOptical,'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off')
    axis off
    latWidth = latLim(2)-latLim(1);
    lonWidth = lonLim(2)-lonLim(1);
    [h2]=plot_2D_energy_slice_geodetic_v1(dataInv, amisrData, magcoords,...
        dataInv(1).energyBin, amisrData.nBeams, timePFISRNo,...
        110, 100,latWidth,lonWidth,false,settings.energySlice.setTimeLabel);

    colormap(axesHandleOptical,'viridis');
    colormap(axesHandleEnergy,'inferno');
    cb1 = colorbar(axesHandleOptical,'eastoutside');
    cb2 = colorbar(axesHandleEnergy,'westoutside');
    ylabel(cb1, '[Rayleigh]');                  
    ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
    caxis(axesHandleOptical,[250 500]);
    caxis(axesHandleEnergy,[8 10]);
    %     caxis(axesHandleEnergy,'auto');
    % Making Energy spectra translucent
    alpha(axesHandleEnergy,0.5);
    % Linking axes together
    set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
    linkaxes([axesHandleOptical,axesHandleEnergy]);

% Storing image
    if setStoreImage == true
        imageName = strcat('figure_',datestr(timeASI(timeIndex),'HH_MM_SS'));
        export_fig(strcat(dirName,'/',storeImageDir,'/',imageName,'.png'),'-r600','-png','-nocrop');
        close(figureHandle);
    end
end