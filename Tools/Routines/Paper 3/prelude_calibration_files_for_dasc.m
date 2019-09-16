%% Download DASC days for calibration

storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir\';
folderStr = 'G:\My Drive\Research\Projects\Paper 3\Data\DASC_Calibration_Files';

timeStr{1} = '14 Dec 2010';
timeStr{2} = '02 Nov 2016';

% download_DASC_FITS(timeStr{1},storeDir);
% fileStr1=find_quite_time_frame(timeStr{1}, storeDir);
%%
% fileStr1 = "G:\My Drive\Research\Projects\Paper 3\Data\WorkDir\20101214\PKR_DASC_0558_20101214_033852.000.FITS";
fileStr1 = "G:\My Drive\Research\Projects\Paper 3\Data\WorkDir\20151112\PKR_DASC_0630_20151112_160443.921.FITS";


% fileStr11 = "G:\My Drive\Research\Projects\Paper 3\Data\WorkDir\20101204\PKR_DASC_0558_20101204_080835.000.FITS";
image1 = fitsread(fileStr1);
%%
% image2 = fitsread(fileStr11);
sensorloc=[65.1260,-147.4789,689 ];
time = fitsfiletimestamp(fileStr1);

%%
AZ_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_Az.FIT'));
EL_2015 = fitsread(strcat(folderStr,filesep,'PKR_DASC_0558_20150213_El.FIT'));

%% Manual Calibration
starCatalogueFitFile = [initialize_root_path,filesep,'energy-height-conversion',...
    filesep,'Tools',filesep,'External Tools',filesep,...
    'skymap',filesep,'hipparcos_extended_catalogue_J2000.fit'];    

stars = get_star_catalogue(starCatalogueFitFile);

[stars.az,stars.el] = RADec2AzEl(rad2deg(stars.RA),rad2deg(stars.DEC),...
    sensorloc(1),sensorloc(2),datestr(time,'yyyy/mm/dd HH:MM:ss'));

realstar = get_actual_stars(stars, 22.5, 4, 0, 0, 0, 1);

%%
toggle = true;

% display_image(toggle,image1,[300 500],'Original Image');

[dasc.az, dasc.el] = get_AzEl_grid(size(image1,1));

%%
indx1 = dasc.el>0;
h1 = figure;
dsign = -1;

    p(1) = plot_DASC_aer(image1(indx1), dasc.az(indx1), dasc.el(indx1), 512, dsign);
    colorbar;
    colormap(get_colormap('k','w'));
    xlim([-120,+120]);
    ylim([-120,+120]);
    hold on;
    
    p(2) = plot_aer_stars(realstar.locationAzEl(:,1), realstar.locationAzEl(:,2),...
        realstar.brightness*80, 'c', 0, 0, 0, dsign);
    
    hold on;
    plot_grid_aer([0, 90], 22.5, 'm');
    caxis([300 500]);
    
    hold on;
    plot_star_names(realstar, 8, 'w', 0, 0, 0, dsign);
    legend([p(1), p(2)], 'Uncalibrated image','Star chart');
    title(datestr(time));
    
%%
    %dx, dy, drot, dsign, dr, k, k0
calPar = [-10, 0, -88, -1, 0, 0, 0]; 
[dasc.azCal, dasc.elCal] = calculate_new_AzEl(dasc.az,dasc.el,calPar);
% dasc.azCal = AZ_2015;
% dasc.elCal = EL_2015;
indx2 = dasc.elCal>0;
dsign = -1;


h2 = figure;
    q(1) = plot_DASC_aer(image1(indx2), dasc.azCal(indx2), dasc.elCal(indx2), 512, dsign);
    colorbar;
    colormap(get_colormap('k','w'));
    xlim([-120,+120]);
    ylim([-120,+120]);
    hold on;
    
    q(2) = plot_aer_stars(realstar.locationAzEl(:,1), realstar.locationAzEl(:,2),...
        realstar.brightness*80, 'c', 0, 0, 0, dsign);
    
    hold on;
    plot_grid_aer([0, 90], 22.5, 'm');
    caxis([300 500]);
    
    hold on;
    plot_star_names(realstar, 8, 'w', 0, 0, 0, dsign);
    legend([q(1), q(2)], 'Calibrated image','Star chart');
    title(datestr(time));
    
    %%
fitswrite(dasc.azCal,...
    "G:\My Drive\Research\Projects\Paper 3\Data\DASC_Calibration_Files\PKR_DASC_20110112_AZ_10deg_Nithin.FITS");
fitswrite(dasc.elCal,...
    "G:\My Drive\Research\Projects\Paper 3\Data\DASC_Calibration_Files\PKR_DASC_20110112_EL_10deg_Nithin.FITS");
% [AZ, EL, accuracy] = calibrate_all_sky_camera(image1, image2, time, sensorloc, [], true);

% download_DASC_FITS(timeStr{2},storeDir);
% fileStr2=find_quite_time_frame(timeStr{1}, storeDir);
% fileStr2

function [newstar, I] = sort_star(newstar)
% Function that sorts stars according to brightness
[newstar.brightness, I] = sort(newstar.brightness,'descend');
newstar.brightness = newstar.brightness';
newstar.location = newstar.location(I,:);
end

function realstar = get_actual_stars(stars,elCutOff,magCutOff,dx,dy,drot, dsign)
if nargin < 7
    dsign = 1;
end

starfilter=stars.vmag<magCutOff & stars.el>elCutOff;
[x,y] = get_aer_stars(stars.az(starfilter), stars.el(starfilter), dx, dy, drot, dsign);
realstar.location = [x, y]; %pixel location
realstar.brightness = stars.relIntensity(starfilter);
realstar.locationAzEl = [stars.az(starfilter), stars.el(starfilter)];
realstar.name = stars.name(starfilter);
[realstar, I] = sort_star(realstar); %s
realstar.locationAzEl = realstar.locationAzEl(I,:);
realstar.name = realstar.name(I);
end


function plot_star_names(sortedStars, n, color, dx, dy, drot ,dsign)

if nargin < 3 || isempty(color)
    color = 'c';
end
if nargin <2 || isempty(n)
    n = 5;
end

hold on;

starNames = string(extractBetween((sortedStars.name(1:1:n)),'(',')'));

for i = 1:1:n
    [x,y] = get_aer_stars(sortedStars.locationAzEl(i,1), sortedStars.locationAzEl(i,2),...
        dx,dy,drot,dsign);
    t=text(x,y,strcat("  ",starNames(i)));
    t.Color = color;
    t.FontSize = 7;
end

end

function p = plot_aer_stars(az,el,relIntensity,colorStr, dx, dy, drot, dsign, dr, k, k0)

if nargin<11
    k0 = 0; %radial distortion
end
if nargin<10
    k = 0; %radial distortion
end
if nargin<9
    dr = 0; % Extent of the radius/ FoV
end
if nargin<8
    dsign = 1;
end

[x,y] = get_aer_stars(az, el, dx, dy, drot, dsign, dr, k, k0);
p = scatter(x,y,relIntensity,colorStr);

end

function [x,y] = get_aer_stars(az, el, dx, dy, drot, dsign, dr, k, k0)
if nargin <9
    k0 = 5; % k0 is an unused parameter.
end
if nargin <8
    k = 0; % k is the distortion parameter
end
if nargin <7
    dr = 0; % Determines the field of view in the image
end
if nargin<6
    dsign = 1;
end
az = rotate_array(az,drot);

r0 = (90-el);
% r1 = radius_with_distortion(el,k,k0);
k = 1+k.*10^-6;
r1 = (90.^((k-1)./k)).*r0.^(1./k);

r = r1 + dr.*r1;
x = dx + dsign.*((r).*sind(az));
y = dy + (r).*cosd(az);
end

function [azNew, elNew] = calculate_new_AzEl(azOld,elOld,x)

[xNew,yNew] = get_aer_stars(azOld,elOld,x(1),x(2),x(3),x(4),x(5),x(6));
azNew = wrapTo360(atan2d(xNew,yNew));
elNew = 90 - xNew./(sind(azNew));
indx = elNew<=0;
azNew(indx) = nan;
elNew(indx) = nan;

end
function [az, el] = get_AzEl_grid(imageSize)

x = linspace(-90,90,imageSize);
y = linspace(-90,90,imageSize);

[X, Y] = meshgrid(x,y);

az = wrapTo360(atan2d(X,Y));
el = 90 - X./(sind(az));

indx = el<=0;
az(indx) = nan;
el(indx) = nan;

end

function f = display_image(toggle,image,clim,titleStr)

if toggle
    
    if nargin<4
        titleStr = '';
    end
    
    f=figure;
    h=pcolor(image);
    set(h,'EdgeColor','none');
    colorbar;
    colormap(get_colormap('k','w'))
    title(titleStr);
    
    if ~(nargin<3) && ~isempty(clim)
        caxis(clim);
    end
end

end


function fileStr = find_quite_time_frame(timeStr, storeDir)

    dirName = strcat(storeDir,datestr(timeStr,'yyyymmdd'));
    [fileStrList] = get_files_in_folder(dirName,'*.FITS');
    fullFileStr = string(strcat(dirName,filesep,fileStrList)');
    
    multiWaitbar('Calculating ...', 0');
    di = 1./length(fullFileStr);
    
    for i = 1:1:length(fullFileStr)
        multiWaitbar('Calculating ...','Increment',di);
        tempData = fitsread(fullFileStr(i));
        intensityArr(i) = nansum(tempData(:));
    end
    multiWaitbar('Calculating ...',1);
    
    [~,I] = min(intensityArr);
    fileStr = fullFileStr(I);
    
end