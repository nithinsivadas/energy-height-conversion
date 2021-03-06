dirName = [initialize_root_path,'LargeFiles\DASC\20110301'];

files = dir(dirName);
fileIndex = find(~[files.isdir]);
fileTempCells  = struct2cell(files);
fileStr = fileTempCells(1,fileIndex);

%% Generating Time Stamp
aldtnum = fitsfiletimestamp(fileStr);
timeASI = (aldtnum-datenum('jan-01-1970'))*(24*3600);
timeASI = unix_to_matlab_time(timeASI);

itimeStart = find_time(timeASI,'01 Mar 2011 9:00:00');
itimeEnd = find_time(timeASI,'01 Mar 2011 10:00:00');

% Creating Image Storage Directory
imageDir ='Figures';  
mkdir(dirName,imageDir);
set(0,'DefaultFigureVisible','off')                                     

    for timeIndex=itimeStart:1:itimeEnd
    ASIDataStr = strcat(dirName,'\',(fileStr(timeIndex)));

    try
        [ASI.dataNew, ASI.lat, ASI.lon, ASI.az_new, ASI.el_new, ASI.sensorloc, ASI.timeDASC] = DASC_aer_to_geodetic...
        (char(ASIDataStr), azOldRes, elOldRes,...
        512, 30, 110);
        catch ME
            warning(['Could not load All sky Camera FITS file: ',fileStr(timeIndex)]);
    end
    %% Create two axes
    figureHandle = figure;
    resize_figure(figureHandle,148,210); %A5 Paper Size
%     display('1');
    % Plot optical data
    [axesHandleOptical, h1]=plot_DASC_geodetic(ASI.dataNew, ASI.timeDASC,...
        ASI.lat, ASI.lon, 512, [63 67], [-153 -143]);
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
       
    [h2]=plot_2D_energy_slice_geodetic_v1(dataInv, amisrData, magcoords,...
        energyBin, nBeams, timePFISRNo, 110, 100,0.4,0.4,false,false);
%     display('2');
    colormap(axesHandleOptical,'viridis');
    colormap(axesHandleEnergy,'inferno');
    cb1 = colorbar(axesHandleOptical,'eastoutside');
    cb2 = colorbar(axesHandleEnergy,'westoutside');
    ylabel(cb1, '[Rayleigh]');                  
    ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
    caxis(axesHandleOptical,[250 500]);
    caxis(axesHandleEnergy,[8 10]);
%     caxis(axesHandleEnergy,'auto');
    alpha(axesHandleEnergy,0.5);
    set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
    linkaxes([axesHandleOptical,axesHandleEnergy]);
    imageName = strcat('figure_',datestr(timeASI(timeIndex),'HH_MM_SS'));
    export_fig(strcat(dirName,'/',imageDir,'/',imageName,'.png'),'-r600','-png','-nocrop');
    close(figureHandle);
    display([num2str(100*(timeIndex-itimeStart)/(itimeEnd-itimeStart)),'%']);
    
    end
% 