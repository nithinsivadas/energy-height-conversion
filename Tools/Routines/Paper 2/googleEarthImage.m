%% Create google earth DASC image
clear all;
%%
h5FileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\Version 2\20080326.001_bc_15sec-full_v3.h5';
figRootFolderStr = 'G:\My Drive\Research\Projects\Paper 2\Submission_v3\Supporting Information\';
figName = 'cover_Art_sphere_2';

%% Load data
dascData = get_2D_plot_inputs_time_independent(h5FileStr,...
    'plotModeStr','OpticalImage');

dascImage = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);

%% Trim data
timeMinStr = '26 Mar 2008 10:30';
timeMaxStr = '26 Mar 2008 11:40';

%%
asi.indx = (find_time(dascData.time,timeMinStr):1:find_time(dascData.time,timeMaxStr))';
asi.data = dascImage(asi.indx,:,:);
asi.lat = permute(reshape(repmat(dascData.latitude(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.lon = permute(reshape(repmat(dascData.longitude(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.az = permute(reshape(repmat(dascData.azimuth(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.el = permute(reshape(repmat(dascData.elevation(:),1,length(asi.indx)),size(asi.data,2),size(asi.data,3),length(asi.indx)),[3,2,1]);
asi.time = repmat(dascData.time(asi.indx),1,size(asi.data,2),size(asi.data,3));

%%
iTimeStr = '26 Mar 2008 11:18';
iTimeIndx = find_time(asi.time(:,1,1),iTimeStr);
x = squeeze(asi.lat(iTimeIndx,:,:));
y = squeeze(asi.lon(iTimeIndx,:,:));
data = squeeze(asi.data(iTimeIndx,:,:));
% interpolate for imagesc
x1 = linspace(62,73,1024);
y1 = linspace(-165,-130,1024);
[X,Y] = meshgrid(x1,y1);
notNan = ~isnan(x);
F = scatteredInterpolant(x(notNan),y(notNan),data(notNan),'linear','none');
data1 = F(X,Y);
nanIndx = isnan(data1);

%%
cLimLow = 340;
cLimHigh = 400;
altitude = 100000;
alphaMatrix = ones(size(data));


kmlFileName = [figRootFolderStr,figName,'.kml'];
imgName = [figRootFolderStr,figName,'.png'];

image = data1;
image(nanIndx) = 0;
% avg = 367.*2^-16;
% n=3;
% sigma = std2(image);
% lowerLim = avg - n*sigma;
% if lowerLim <=0
%   lowerLim=0;
% end
% upperLim = avg + n*sigma;
% % if upperLim >=1
%   upperLim=1;
% end
% image = imadjust(image, [lowerLim upperLim],[]);
alphaMatrix(flipud(image')<350) = 0;

figure;
imagesc(y1,x1,image',[cLimLow,cLimHigh]);

% cmap = get_colormap('k','g');
cmap = colormap(viridis);
colorbar;

output = ge_imagesc(y1,x1,flipud(image'),...
                    'imgURL',imgName,...
                   'cLimLow',cLimLow,...
                  'cLimHigh',cLimHigh,...
                  'altitude',altitude,...
              'altitudeMode','absolute',...
                  'colorMap',cmap,...
               'alphaMatrix',alphaMatrix);

ge_output(kmlFileName,output,'name',kmlFileName);