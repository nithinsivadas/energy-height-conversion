%% Configuring the paths

if strcmp(get_computer_name,'nithin-carbon')
    dataDir = 'G:\My Drive\Research\Projects\Data\';
    storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\';
elseif strcmp(get_computer_name,'aurora1-optiplex-780')
    dataDir = '/media/nithin/PFISR_002_006/Nithin/Data/';
    storeDir = '/media/nithin/PFISR_002_006/Nithin/Data/Paper_3/';
elseif strcmp(get_computer_name,'scc-lite')
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/scratch/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
else
    dataDir = '/projectnb/semetergrp/nithin/Data/';
    storeDir = '/projectnb/semetergrp/nithin/Data/Paper_3/';
    jobDir = '/projectnb/semetergrp/nithin/local_cluster_object';
    if exist(jobDir,'dir')==7
        system(['rm -rf ',jobDir,filesep,'*']); 
    end
end

%% Getting input files
workDir = strcat(storeDir,'WorkDir');
fileNameList = struct2cell(dir([workDir,strcat(filesep,'*_dascData.h5')]));
filePathStr = strcat(strcat(workDir,filesep),string(fileNameList(1,:)'));
h5Calibration = fullfile(storeDir,'DASC_Calibration_Files','dasc_calibration.h5');
if isempty(filePathStr)
    error('No files in WorkDir');
end
if ~isfile(h5Calibration)
    error(['Calibration file not available: ',h5Calibration]);
end
calibration = read_h5_data(h5Calibration);
%%
setSample = true; %% Plot only samples - 99 frames for each substorm

%%
myCluster = parcluster("local");
myCluster.JobStorageLocation = jobDir;
try
    delete(myCluster.Jobs);
catch
end

for i=1:1:length(filePathStr)
    j(i)=batch(myCluster, @batch_process,0,{i,filePathStr,calibration,workDir,setSample});
%     batch_process(i,filePathStr,calibration,workDir,setSample);
end

%% Wait for each job to finish before quitting matlab
for i=1:1:length(filePathStr)
    wait(j(i));
end

function batch_process(i,filePathStr, calibration, workDir, setSample)
    fileName = filePathStr(i);
    tempStr = strsplit(fileName,filesep);
    tempStr1 = strsplit(tempStr(end),'.');
    videoFileName = strcat(tempStr1(1),'.avi');
    imageDir = ['images_',num2str(i)];
    imageLongDir = fullfile(tempStr{1:end-1},imageDir);
    if isunix
        imageLongDir = char(strcat(filesep,imageLongDir));
    end
    [status,~] = create_images(fileName,imageDir, true, calibration, setSample); 
    create_video(workDir,imageDir,videoFileName);
    try
    rmdir(imageLongDir,'s');
    catch
    end
end


function [status,imageDir] = create_images(fileName, imageDirName, setStoreImage,...
    calibration, setSample)

status = 'Failed';
tempStr = strsplit(fileName,filesep);
tempStr1 = strsplit(tempStr{end},'_');
pfisrSearch = fullfile(tempStr{1:end-1},strcat(tempStr1{1},"*_pfisrData.h5"));
if isunix
    pfisrSearch = char(strcat(filesep,pfisrSearch));
end
pfisrFile = struct2cell(dir(pfisrSearch));
if ~isempty(pfisrFile(1,:))
   if length(pfisrFile(1,:))>1
       warning('More than one PFISR data file, taking just the first one.');
   end
   h5pfisr = fullfile(pfisrFile{2,1},pfisrFile{1,1});
   pfisrData.time = h5read(h5pfisr,'/time');
   pfisrData.energyBin = h5read(h5pfisr,'/energy/energyBin');
   pfisrData.energyFlux = h5read(h5pfisr,'/energy/energyFlux');
   pfisrData.alt = h5read(h5pfisr,'/alt');
   pfisrData.az = h5read(h5pfisr,'/az');
   pfisrData.el = h5read(h5pfisr,'/el');
   pfisrData.electronDensity = h5read(h5pfisr,'/inputData/Ne');
   pfisrData.sigmaH = h5read(h5pfisr,'/conductivity/hall');
   pfisrData.sigmaP = h5read(h5pfisr,'/conductivity/pederson');
   pfisrData.lat = h5read(h5pfisr,'/lat');
   pfisrData.lon = h5read(h5pfisr,'/lon');
   pfisrData.point.x = (90-pfisrData.el(1)).*sind(pfisrData.az(1));
   pfisrData.point.y = (90-pfisrData.el(1)).*cosd(pfisrData.az(1));
else
   pfisrData = [];
end
    
[az,el,lat,lon] = get_calibration_coordinates(calibration,tempStr1{2});
[Fae,~,flagPixel] = get_scattered_points(az,el,lat,lon);

imageDir = strcat(fullfile(tempStr{1:end-1}),filesep,imageDirName);
if isunix
    imageDir = char(strcat(filesep,imageDir));
end

data = read_h5_data(fileName);
indx=find(strcmp(string(data.Name),'ASI'));

    if(~isfolder(imageDir))
        mkdir(imageDir);
    end
    

    for i=1:1:length(indx)
        tempStr = strsplit(data.Path{indx(i)},'/');
        wavelengthStr(i) = tempStr(2);
        ASI{i} = data.Data{indx(i)};
        time{i} = unix_to_matlab_time(data.Data{strcmp(string(data.Path),strcat('/',tempStr(2),'/','time'))});
        tLength(i) = length(time{i});
        tempArr = ASI{i}(:);
        tempArr(tempArr>10000) = nan;
        asiPlotVar(i).median = nanmedian(tempArr);
        asiPlotVar(i).mean = nanmean(tempArr);
        asiPlotVar(i).std = nanstd(tempArr);
    end

[~,maxTimeIndx] = max(tLength);
mainTime = time{maxTimeIndx};
    if setSample
        timeArr = round(linspace(1,length(mainTime),min(99,length(mainTime)))); %only 99 frames per sample
    else
        timeArr = 1:1:length(mainTime);
    end

    for iTime = timeArr
        h=figure('visible','off');
        h=create_figure(h,Fae,flagPixel,datestr(mainTime(iTime)),...
            time,ASI,asiPlotVar,wavelengthStr,pfisrData);         
        if setStoreImage == true
           export_fig(fullfile(imageDir,strcat('Figure_',datestr(mainTime(iTime),'HH_MM_ss'),'.png')),...
               '-r300','-png','-nocrop');
           close(h);
        end
    end
    
end

function h=create_figure(h,...
    Fae,flagPixel,timeStr,time,ASI,asiPlotVar,wavelengthStr,...
    pfisrData)
  
    p = panel();

    nLambda= length(wavelengthStr);

    if nLambda == 3
        nX = 2;
        nY = 2;
        resize_figure(h,210,250);
        panelHeight = 100; %mm
        panelBreadth = 100; %mm
        rowHeights = repmat({{panelHeight}},1,nY);
        colBreadths = repmat({{panelBreadth}},1,nX);
        p.pack(rowHeights, colBreadths);
        p.de.margin = 4;
               
        if ~isempty(pfisrData)
            p(1,2).marginleft=30;
            try
            p=plot_pfisr_panels(p, pfisrData, timeStr);
            catch
                p(1,2).select();
                set(gca,'XColor','none','YColor','none');
            end
        else
            p(1,2).select();
            set(gca,'XColor','none','YColor','none');
        end
        
        iLambda = (strcmp(wavelengthStr,'0558'));
        p(1,1).select();
        plot_dasc_color(time,timeStr,Fae,flagPixel,ASI,...
            iLambda,wavelengthStr,asiPlotVar,pfisrData,'y',1);     
        
        iLambda = (strcmp(wavelengthStr,'0428'));
        p(2,1).select();
        plot_dasc_color(time,timeStr,Fae,flagPixel,ASI,...
            iLambda,wavelengthStr,asiPlotVar,pfisrData,'c',2);
        
        iLambda = (strcmp(wavelengthStr,'0630'));
        p(2,2).marginleft=15;
        p(2,2).select();
        plot_dasc_color(time,timeStr,Fae,flagPixel,ASI,...
            iLambda,wavelengthStr,asiPlotVar,pfisrData,'m',3);
        
    elseif nLambda == 1
        
        nX = 2;
        nY = 1;
        panelHeight = 100; %mm
        panelBreadth = 100; %mm
        rowHeights = repmat({{panelHeight}},1,nY);
        colBreadths = repmat({{panelBreadth}},1,nX);
        p.pack(rowHeights, colBreadths);
        p.de.margin = 4;
                 
        if ~isempty(pfisrData)
            resize_figure(h,110,250);
            p(1,2).marginleft=30;
            p=plot_pfisr_panels(p, pfisrData, timeStr);
        else
            p(1,2).select();
            set(gca,'XColor','none','YColor','none');
            resize_figure(h,110,135);
        end
        
        iLambda = 1;
        p(1,1).select();
        plot_dasc_color(time,timeStr,Fae,flagPixel,ASI,...
            iLambda,wavelengthStr,asiPlotVar,pfisrData,'y',1);
    end
    
end


function plot_dasc_color(time,timeStr,Fae,flagPixel,ASI,...
    iLambda,wavelengthStr,asiPlotVar,pfisrData,colorStr,turn)
        
        set(gca,'XColor','none','YColor','none');
        try
        closestTimeIndx = find_time(time{iLambda},timeStr);
        image = ASI{iLambda}(closestTimeIndx,:,:);
        imean = nanmean(image(:));
        istd = nanstd(image(:));
        
        
        plot_DASC_aer_v1(Fae, flagPixel, image);
        if turn==1
        text(0.5,1,{datestr(time{iLambda}(closestTimeIndx),'dd mmm yyyy')},...
            'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom');
        end
        text(0.9,0,{datestr(time{iLambda}(closestTimeIndx),'HH:MM:ss'),strcat(wavelengthStr(iLambda),' nm')},...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom');
        colormap(gca,get_colormap('k',colorStr));
        
        cmax = min(ceil(asiPlotVar(iLambda).mean+asiPlotVar(iLambda).std),imean+4*istd);
        
        caxis([floor(asiPlotVar(iLambda).median), cmax]);       
        c = colorbar_thin();
        set(c, 'YAxisLocation', 'left');
        if ~isempty(pfisrData)
            hold on; plot(pfisrData.point.x,pfisrData.point.y,'ow');
        end
        catch             
        end
end
function p=plot_pfisr_panels(p, pfisrData, timeStr)
            
            nX = 1;
            nY = 4;
            panelHeight = 20; %mm
            panelBreadth = 80; %mm
            rowHeights = repmat({{panelHeight}},1,nY);
            colBreadths = repmat({{panelBreadth}},1,nX);
            p(1,2).pack(rowHeights, colBreadths);
            p(1,2).de.margin = 4;
            
            closestTimeIndxPFISR = find_time(pfisrData.time,timeStr);
            x = [pfisrData.time(closestTimeIndxPFISR) , pfisrData.time(closestTimeIndxPFISR)];
            
            p(1,2,1,1).select();
            plot_2D_time_series(pfisrData.time,pfisrData.alt,pfisrData.electronDensity,0.5);
            c11=colorbar_thin();
            colormap(gca,magma);
            ylabel(c11,'[m^-^3]');
            set(gca,'colorscale','log');
            ylabel({'N_e','[km]'});
            yAlt = [60, 140];
            ylim(yAlt);
            caxis([10^8, 10^12]);
            c11.Ticks=[10^8,10^10,10^12];
            hold on;
            line(x,yAlt,'Color','cyan');
            title('PFISR');
            
            p(1,2,2,1).select();
            plot_2D_time_series(pfisrData.time,pfisrData.energyBin,pfisrData.energyFlux,0.5);
            c21=colorbar_thin();
            colormap(gca,inferno);
            ylabel(c21,'[ev/ m^2 sr s eV]');
            set(gca,'colorscale','log','YScale','log',...
                'YTick',[10, 30, 100, 500]*1000,'YTickLabel',{'10','30','100','500'});
            ylabel({'\phi(E)','[keV]'});
            yEnergy = [5, 900].*1000;
            ylim(yEnergy);
            caxis([10^8, 10^12]);
            c21.Ticks=[10^8,10^10,10^12];
            hold on;
            line(x,yEnergy,'Color','cyan');
            
            p(1,2,3,1).select();
            plot_2D_time_series(pfisrData.time,pfisrData.alt,pfisrData.sigmaH,0.5);
            c31=colorbar_thin();
            colormap(gca,viridis);
            ylabel(c31,'[S/m]');
            set(gca,'colorscale','log');
            ylabel({'\sigma_H','[km]'});
            ylim(yAlt);
            caxis([10^-6, 10^-2]);
            c31.Ticks=[10^-6,10^-4,10^-2];
            hold on;
            line(x,yAlt,'Color','m');
            
            p(1,2,4,1).select();
            plot_2D_time_series(pfisrData.time,pfisrData.alt,pfisrData.sigmaP,0.5);
            c41=colorbar_thin();
            colormap(gca,viridis);
            ylabel(c41,'[S/m]');
            set(gca,'colorscale','log');
            ylabel({'\sigma_P','[km]'});
            ylim(yAlt);
            caxis([10^-7, 10^-3]);
            c41.Ticks=[10^-7,10^-5,10^-3];
            label_time_axis(pfisrData.time,true,0.5);
            hold on;
            line(x,yAlt,'Color','m');  

end

function h2=plot_DASC_aer_v1(Fae, flagPixel, dataNew, imageSize)

if nargin<4
    %imageSize = max(size(dataNew));
    imageSize = 512; %%?? Only as long as we don't have the right calibration files
    if max(size(dataNew))>imageSize
        dataNew = modify_matrix_size(squeeze(dataNew),imageSize,imageSize);
    end
end

xq = linspace(-90,90,imageSize);
yq = linspace(-90,90,imageSize);
[X,Y] =ndgrid(xq,yq);

Fae.Values = dataNew(flagPixel);
Vq1 = Fae(X, Y);

h2=pcolor(X,Y,Vq1);
set(h2,'EdgeColor','none');
% colorbar;
axis equal
% hold on; annotation('textarrow',[0.15 0.20],[0.8 0.8],'string','W');
% hold on; annotation('textarrow',[0.2 0.2],[0.85 0.90],'string','N');
hold on; text(-77,75,'N');
hold on; line([-75 -75],[70 60],'color','r');
hold on; text(-77,55,'S');

hold on; text(-69,65,'W');
hold on; line([-80 -70],[65 65],'color','r');
hold on; text(-85,65,'E');

end

function [Fae,Fll,flagPixel] = get_scattered_points(az,el,lat,lon,elMin)

if nargin<5 || isempty(elMin)
    elMin = 15;
end

nanIndx = isnan(lat)|isnan(lon);
flagPixel = el>elMin & ~nanIndx;

x = (90-el).*sind(az);
y = (90-el).*cosd(az);

Fae = scatteredInterpolant(x(flagPixel),y(flagPixel),zeros(length(x(flagPixel)),1),'linear','none');
Fll = scatteredInterpolant(lat(flagPixel),lon(flagPixel),zeros(length(x(flagPixel)),1),'linear','none');

end


function [az,el,lat,lon,minElFlag] = get_calibration_coordinates(calibration, timeStr)
     
    timeRange = calibration.Data{strcmp(string(calibration.Name),'timeRange')};
    timeSample = datenum(timeStr,'YYYYmmdd');
    
    chosenFileIndx = timeRange(timeRange(:,2)<=timeSample & timeRange(:,3)>=timeSample,1);
    
    azDataStr = strcat('/',num2str(chosenFileIndx),'/AZ');
    elDataStr = strcat('/',num2str(chosenFileIndx),'/EL');
    latDataStr = strcat('/',num2str(chosenFileIndx),'/lat');
    lonDataStr = strcat('/',num2str(chosenFileIndx),'/lon');
    minElFlagStr = strcat('/',num2str(chosenFileIndx),'/minElFlag');
    
    az = calibration.Data{strcmp(string(calibration.Path),azDataStr)};
    el = calibration.Data{strcmp(string(calibration.Path),elDataStr)};
    lat = calibration.Data{strcmp(string(calibration.Path),latDataStr)};
    lon = calibration.Data{strcmp(string(calibration.Path),lonDataStr)};
    minElFlag = calibration.Data{strcmp(string(calibration.Path),minElFlagStr)};
    
    if chosenFileIndx==1
        if chosenFileIndx == 2
            az = 180-az';
        else
            az = az';
        end
        el = el';
        lat = lat';
        lon = lon';
        minElFlag = minElFlag';
    end
end
