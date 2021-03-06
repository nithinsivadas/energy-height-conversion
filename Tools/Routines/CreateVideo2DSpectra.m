clear all;
close all;

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Data/2D_energy_spectra_26032008_800_to_1300Hr.mat');
azStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_AZ_10deg.FITS';
elStr='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26/PKR_Cal_before_2011/PKR_20111006_EL_10deg.FITS';
include_new_colormaps;

% Creating DASC file list and time array
dirName='/home/nithin/Documents/git-repos/Largefiles/PokerFlat_DASC_08_03_26';
files = dir(dirName);
fileIndex = find(~[files.isdir]);
fileTempCells  = struct2cell(files);
fileStr = fileTempCells(1,fileIndex);

%% Generating Time Stamp
aldtnum = fitsfiletimestamp(fileStr);
timeASI = (aldtnum-datenum('jan-01-1970'))*(24*3600);
timeASI = unix_to_matlab_time(timeASI);

itimeStart = find_time(timeASI,'26 Mar 2008 11:00:00');
itimeEnd = find_time(timeASI,'26 Mar 2008 12:00');

% Creating Image Storage Directory
imageDir ='Figures';  
mkdir(dirName,imageDir);
set(0,'DefaultFigureVisible','off')                                     

for timeIndex=itimeStart:1:itimeEnd
    ASIDataStr = strcat(dirName,'/',(fileStr(timeIndex)));


    azOldRes=fitsread(azStr);
    elOldRes=fitsread(elStr);
    try
        [dataNew, lat, lon, az_new, el_new, sensorloc, timeDASC] = DASC_aer_to_geodetic...
        (char(ASIDataStr), azOldRes, elOldRes,...
        512, 30, 110);
    catch 
        continue;
    end
    %% Create two axes
    figureHandle = figure;
    resize_figure(figureHandle,148,210); %A5 Paper Size
    
    % Plot optical data
    [axesHandleOptical, h1]=plot_DASC_geodetic(dataNew, timeDASC, lat, lon, 512, [63 67], [-152 -143]);
    
    % Identify the closest time in the PFISR data to the current optical
    % time
    timePFISRNo = find_time(data(1).time,datestr(timeASI(timeIndex)));
    % Plot energy data
    axesHandleEnergy = axes;
    axesm('lambertstd','MapLatLimit',getm(axesHandleOptical,'MapLatLimit'),...
            'MapLonLimit',getm(axesHandleOptical,'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off')
    axis off
       
    [h2]=plot_2D_energy_slice(data, magcoords, energyBin, nBeams, timePFISRNo, 110, 100, false);

    colormap(axesHandleOptical,'viridis');
    colormap(axesHandleEnergy,'inferno');
    cb1 = colorbar(axesHandleOptical,'eastoutside');
    cb2 = colorbar(axesHandleEnergy,'westoutside');
    ylabel(cb1, '[Rayleigh]');                  
    ylabel(cb2, 'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
    caxis(axesHandleOptical,[0 300]);
    caxis(axesHandleEnergy,[8 10]);
%     caxis(axesHandleEnergy,'auto');
    alpha(axesHandleEnergy,0.5);
    set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'));
    linkaxes([axesHandleOptical,axesHandleEnergy]);
    imageName = strcat('figure_',datestr(timeASI(timeIndex),'HH_MM_SS'));
    export_fig(strcat(dirName,'/',imageDir,'/',imageName,'.png'),'-r600','-png','-nocrop');
    close(figureHandle);
    display([num2str(100*(timeIndex-itimeStart)/(itimeEnd-itimeStart)),'%']);
    
end;
% 
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample Images/Sample_aurora_energy.pdf' -pdf -nocrop
% save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample_aurora_energy.svg');
% export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 2/Sample Images/Sample_aurora_energy.png' -r600 -png -nocrop
% 

%% Create Video
display('Generating video file');
outputVideo = VideoWriter(fullfile('/home/nithin/Documents/git-repos/Largefiles/26_Mar_2008_Optical_and_Energy_01.avi'));
outputVideo.FrameRate = 8;
open(outputVideo)                                                                                                                                        
imageFiles = dir(strcat(dirName,'/',imageDir));
imageFileIndex = find(~[imageFiles.isdir]);
fileTempCells  = struct2cell(imageFiles);
imageFileStr = fileTempCells(1,imageFileIndex);
for ii = 1:length(imageFileStr)
   img = imread(fullfile(dirName,imageDir,imageFileStr{ii}));
   writeVideo(outputVideo,img)
end                                                                                                                                    
close(outputVideo);
set(0,'DefaultFigureVisible','on')       
