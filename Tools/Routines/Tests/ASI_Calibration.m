%% Calibrate Poker Flat ASC
clear all;

stars=get_star_catalogue;
% fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';
fileStr = 'C:\Users\nithin\Downloads\20080326.001_bc_15sec-full_v2.h5';
%% Step 1: Estimate Darkest hour
[totalIntensity,timeArr]=estimate_darkest_frame(fileStr);
[val,indx]=min(totalIntensity);
timeStr = datestr(timeArr(indx));
time = timeArr(indx);

dasc.sensorLoc = h5read(fileStr,'/DASC/sensorloc');

polarisIndx = find(stars.HIP==11767); % Validation star : Polaris

% Approximation
[stars.az,stars.el] = RADec2AzEl(rad2deg(stars.RA),rad2deg(stars.DEC),...
    dasc.sensorLoc(1),dasc.sensorLoc(2),datestr(time,'yyyy/mm/dd HH:MM:ss'));

% Sophesticated calculation of Polaris
[polaris.azAccurate,polaris.elAccurate] = get_star_az_el...
    (stars.RA(polarisIndx),stars.DEC(polarisIndx),...
    stars.pmRA(polarisIndx),stars.pmDEC(polarisIndx),stars.parallax(polarisIndx),...
    stars.RV(polarisIndx),time,deg2rad(dasc.sensorLoc(1)),deg2rad(dasc.sensorLoc(2)),dasc.sensorLoc(3));

% From Stellaris - cross-check
polaris.azStellaris = 359+17./60+39.1./36000;
polaris.elStellaris = 64+30./60+7.8./3600;

% From Approximation
polaris.az = stars.az(polarisIndx);
polaris.el = stars.el(polarisIndx);

%% Reading DASC data
dasc.ASI = read_h5_variable_at_time_v2(fileStr,'/DASC/ASI',[],timeStr);
dasc.az = (h5read(fileStr,'/DASC/azCalData'))';
dasc.el = (h5read(fileStr,'/DASC/elCalData'))';
%% Step 2: Identify Hot Pixles
% Load fits File
asiPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC\20080326\';
asi1File = 'PKR_DASC_0000_20080326_103958.000.FITS';
asi2File = 'PKR_DASC_0000_20080326_104018.000.FITS';
asi3File = 'PKR_DASC_0000_20080326_133018.000.FITS';

image1 = fitsread([asiPath,asi1File]);
ASI2 = fitsread([asiPath,asi2File]);
image2 = fitsread([asiPath,asi3File]);


[hotPixels] = identify_hot_pixels(image1, image2, 10);

%% Step 3: Remove background
[image1BkgRemRow, backgroundRow, imRowNoise, pRow, muRow, a] = remove_background(image1,5);
[image1BkgRem,background1, imColNoise] = remove_background(image1BkgRemRow');
image1BkgRem = image1BkgRem';
imColNoise = imColNoise';
background1 = background1';
%%
figure; 
imagesc(image1BkgRemRow);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed (Row fit)');
figure; 
imagesc(imRowNoise);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 10]); title('ASI noise (Row fit)');
figure;
imagesc(background1);
colorbar;
colormap(get_colormap('k','w'));
caxis([200 800]); title('Background with stars removed (Row fits)');
%%
figure; 
imagesc(image1BkgRem);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed (Col fit)');
%%
figure;
imagesc(imColNoise);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 10]); title('ASI noise (Col fit)');
figure;
imagesc(background1);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 50]); title('Background with stars removed (Col fits)');

%% Step 4: Calculate the sigma_n / noise level
fnoise = @(x) nanstd(x(:));
totalNoise = imColNoise+imRowNoise;
totalNoise(totalNoise==0)=nan;
sigma_n = nlfilter(totalNoise,[9 9],fnoise);
%%
figure; 
imagesc(sigma_n);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 10]); title('\sigma_n Noise at each pixel');

%% Step 5: Noise spike removal
image1NoiseRem = remove_noise_spikes(image1BkgRem, sigma_n);
figure; 
imagesc(image1NoiseRem);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed + noise spikes removed (Col fit)');

%% Step 6: Star extraction 
[starImage] = faint_star_extracter(image1BkgRem, sigma_n);
% Rough removal of noise below 20 pixels
%
starImage(starImage <= 20) = 0;
starImage(1:125,:) = 0;
starImage(1000:1024,:) = 0;
%%
% faintStarImageModified = modify_matrix_size(starImage, 512, 512);
%%
[brightStarImage] = bright_star_extracter(ASINew1, sigma_n);
% Rough removal of noise below 20 pixels
%
brightStarImage(brightStarImage <= 20) = 0;
brightStarImage(1:125,:) = 0;
brightStarImage(1000:1024,:) = 0;

%% Step 7: Extracting the star position and brightness

imstarStruct = extract_stars(starImage);

% Process star position and brightness for astrometry
astroFileStr = 'G:\My Drive\Research\Projects\Paper 3\Data\StarCalibration\sivadas_astrometry_01.txt';
dascstar = astrometry_file(imstarStruct,22.5);

%%
figure; 
h=pcolor(starImage);
set(h,'EdgeColor','none');
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('Faint stars (extraction algorithm)');
hold on;
scatter(dascstar.location(:,1),dascstar.location(:,2),50*dascstar.brightness/max(dascstar.brightness),'r');
%%
figure; 
imagesc(brightStarImage);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('Bright stars (extraction algorithm)');
%%
i = 50;
x=a(i).fitRangeIndx;
figure;
plot(x,polyval(pRow(i,:),x));
hold on;
plot(x,image1(i,x));
% ylim([0,1000]);
%%
% change scale of existing az & el
modAz = modify_matrix_size(dasc.az,1024,1024);
modEl = modify_matrix_size(dasc.el,1024,1024);
%%
realstar = get_actual_stars(stars,22.5, 4, 0, 0, 0, 1); %some problem with realstar
%% Plotting
figure;
colormap(viridis);
% Plotting all-sky camera image
indx = modEl>0;
plot_DASC_aer(image1(indx), rotate_array(modAz(indx),0), modEl(indx), 1024, -1);
hold on; 
plot_grid_aer([0,90],[22.5],'m');
% colorbar;
xlim([-120,+120]);
ylim([-120,+120]);
caxis([300 450]);
hold on;

% Stars
dx = 8;
dy = +2.5;
drot = +90;

% starfilter=stars.vmag<4 & stars.el>22.5;
% [x,y] = get_aer_stars(stars.az(starfilter),stars.el(starfilter), dx, dy, drot, -1);
% scatter(x,y,stars.relIntensity(starfilter)*50,'r');
[x,y] = get_aer_stars(realstar.locationAzEl(:,1), realstar.locationAzEl(:,2), dx, dy, drot, -1);
scatter(x,y,realstar.brightness*30,'r');
hold on;

plot_aer_stars(polaris.azStellaris,polaris.elStellaris,10,'w', dx, dy, drot, -1);
hold on;
plot_aer_star_label(polaris.az, polaris.el, 'Polaris', [],  dx, dy, drot, -1);
hold on;
plot_aer_stars(polaris.azAccurate, polaris.elAccurate,10,'c', dx, dy, drot, -1);

title(timeStr);
%%


%%
tic
modAz1=modAz;
modAz1(modAz1==0) = nan;
modEl1=modEl;
modEl1(modEl1==0) = nan;

[dascstar1,x,fval] = calibrate_stars(realstar,dascstar,modAz1, modEl1);
toc
x
fval
%
%%
figure; 
indx = modEl1>0;
plot_DASC_aer(starImage(indx), modAz1(indx), modEl1(indx), 1024, 1);
colorbar;
colormap(get_colormap('k','w'));
caxis([0,100]);
hold on;
plot_aer_stars(dascstar1.locationAzEl(:,1),dascstar1.locationAzEl(:,2), dascstar1.brightness*20,'r',0,0,0,1);

%%
figure;
% indx = modEl>0;
% plot_DASC_aer(image1(indx), rotate_array(modAz(indx),0), modEl(indx), 1024, 1);
% hold on;
plot_aer_stars(realstar.locationAzEl(:,1),realstar.locationAzEl(:,2),realstar.brightness*50,'r',0,0,0,1);
hold on;
plot_aer_stars(dascstar1.locationAzEl(:,1),dascstar1.locationAzEl(:,2),dascstar1.brightness*10,'b',x(1),x(2),x(3), x(4));



[azdasc1, eldasc1] = calculate_new_AzEl(dascstar1.locationAzEl(:,1),dascstar1.locationAzEl(:,2),x);
hold on;
plot_aer_stars(azdasc1, eldasc1, dascstar1.brightness*10,'g',0,0,0,1);
%% Everything is working now; 
%%
figure;
colormap(viridis);
% Plotting all-sky camera image
[azNew, elNew] = calculate_new_AzEl(modAz1,modEl1,x);
indx = elNew>0;
plot_DASC_aer(image1(indx), azNew(indx),elNew(indx), 1024, -1);
hold on; 
plot_grid_aer([0,90],[22.5],'m');
% colorbar;
xlim([-120,+120]);
ylim([-120,+120]);
caxis([300 450]);
hold on;

plot_aer_stars(azdasc1,eldasc1,dascstar1.brightness*50,'r',0,0,0,-1);

hold on;
plot_aer_stars(realstar.locationAzEl(:,1),realstar.locationAzEl(:,2),realstar.brightness*50,'c',0,0,0,-1);
title(timeStr);

%% Please carefully replicate it


function [azNew, elNew] = calculate_new_AzEl(azOld,elOld,x)
%     azOld(azOld==0)=nan;
%     elOld(elOld==0)=nan;
    [xNew,yNew] = get_aer_stars(azOld,elOld,x(1),x(2),x(3),x(4));
%     problem.objective = @solveAzEl;
%     problem.x0 = [0,90];
%     problem.solver = 'fsolve';
%     problem.options = optimoptions('fsolve','Display','none');
%     multiWaitbar('Solving Az El Values',0);
    azNew = wrapTo360(atan2d(xNew,yNew));
    elNew = 90 - xNew./(sind(azNew));
    
%     di = 1./length(xNew);
%     for i = 1:1:length(xNew)
%         for j = 1:1:length(yNew)
%             xN = xNew(i,j);
%             yN = yNew(i,j);
%             if ~isnan(xN) && ~isnan(yN)
%                 a = fsolve(problem);
%                 azNew(i,j) = a(1);
%                 elNew(i,j) = a(2);
%             else
%                 azNew(i,j) = nan;
%                 elNew(i,j) = nan;
%             end
%         end
%         multiWaitbar('Solving Az El Values','Increment',di);
%         disp(i);
%     end
%     
%     function F = solveAzEl(a)
%         a(1) = rotate_array(a(1),x(3));
%         F(1) = x(4).*(90-a(2)).*sind(a(1)) - xN + x(1);
%         F(2) = (90-a(2)).*cosd(a(1)) - yN + x(2);
%     end
%     
%     azNew = [];
%     elNew = [];
end

function [dascstar,x,fval] = calibrate_stars(realstar,dascstar,azOld,elOld)
    
    dloc = round(dascstar.location);
    dascstar.brightness = dascstar.brightness./max(dascstar.brightness);
    lindx = sub2ind(size(azOld),dloc(:,2),dloc(:,1)); %% There was issue, x - column , y - are rows
    dascstar.locationAzEl = [azOld(lindx), elOld(lindx)];
    minElFilter = find(dascstar.locationAzEl(:,2)>=22.5);
    dascstar.locationAzEl = dascstar.locationAzEl(minElFilter,:);
    dascstar.brightness = dascstar.brightness(minElFilter);
    dascstar.location = dascstar.location(minElFilter,:);
    ndasc = length(dascstar.brightness);
    
    x0 = [0,0,-90,1];
    nvars = 4;
    lb = [-10,-10,-180, -1];
    ub = [+10,+10,+180,+1];
    IntCon = 4;
    
    pflag = 0;
    nflag = 0;
    k = 0;
    while(pflag == 0 || nflag ==0)
        k = k+1;
        [y(k,:),fval(k),exitflag] = ga(@starDistance,nvars,[],[],[],[],lb,ub,[],IntCon);
        
        if y(k,4) == -1
            nflag = 1;
        end
        
        if y(k,4) == 1
            pflag = 1;
        end
    end
    
    [~,minIndx] = min(fval);
    x = y(minIndx,:);    
    
    function dmin=starDistance(x)
        [x1,y1] = get_aer_stars(dascstar.locationAzEl(:,1),dascstar.locationAzEl(:,2),x(1),x(2),x(3),x(4));
        dascstar.newxy = [x1,y1];
        D = pdist2(dascstar.newxy,realstar.location(1:ndasc,:));
        dmin = sum(min(D));
    end

end
% 
% function [x1, y1, az1, el1] = transform_aer(az, el, del, drot, dsign)
% 
%     az1 = rotate_array(az, drot);
%     el1 = el1-
% 
% end


function realstar = get_actual_stars(stars,elCutOff,magCutOff,dx,dy,drot, dsign)
    if nargin < 7
        dsign = 1;
    end

%     stars.el = abs(stars.el);
    starfilter=stars.vmag<magCutOff & stars.el>elCutOff;
    [x,y] = get_aer_stars(stars.az(starfilter), stars.el(starfilter), dx, dy, drot, dsign);
    realstar.location = [x, y];
    realstar.brightness = stars.relIntensity(starfilter);
    realstar.locationAzEl = [stars.az(starfilter), stars.el(starfilter)];
    [realstar, I] = sort_star(realstar);
    realstar.locationAzEl = realstar.locationAzEl(I,:);
end

function [x,y] = get_aer_stars(az,el,dx,dy,drot,dsign)
    if nargin<6
        dsign = 1;
    end
    az = rotate_array(az,drot);
    x = dx+dsign*(90-el).*sind(az);
    y = dy+(90-el).*cosd(az);  
end

function plot_aer_stars(az,el,relIntensity,colorStr, dx, dy, drot, dsign)
    if nargin<8
        dsign = 1;
    end
    [x,y] = get_aer_stars(az, el, dx, dy, drot, dsign);
    scatter(x,y,relIntensity,colorStr);
end

function plot_aer_star_label(az,el,textString,colorStr, dx, dy, drot, dsign)
    if nargin<8
        dsign = 1;
    end

    [x,y] = get_aer_stars(az, el, dx, dy, drot, dsign);
    
    if nargin<4 || isempty(colorStr)
        colorStr = 'r';
    end
    
    text(x, y, textString, 'Color', colorStr);
    
end

function [totalIntensity,timeArr]=estimate_darkest_frame(h5FileStr)

asi = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
timeArr = unixtime2matlab((h5read(h5FileStr,'/DASC/time'))');
totalIntensity = sum(sum(asi,3),2);

end

function [hotPixels] = identify_hot_pixels(image1, ASI2, threshold)
    
    if nargin<3
        threshold = 2;
    end
    
    temp = image1 - ASI2;
    temp2 = temp;
    hotPixels=zeros(size(temp2));
    hotPixels(abs(temp2)<=threshold) = 1;
    
end

function [ASINew, background, ASINoise, pRow, muRow, a] = remove_background(ASI, nPoly)

% Remove background along the row

if nargin < 2
    nPoly = 3;
end

% ASINew, background
rowIndx = 1:1:size(ASI,2);
colIndx = 1:1:size(ASI,1);
ASINew = zeros(size(ASI));
background = zeros(size(ASI));
ASItemp = ASI;
ASINoise = zeros(size(ASI));
pRow = zeros([length(rowIndx),nPoly+1]);
% First iteration
for ifit = 1:2
    for i = rowIndx
        lchord = sqrt(512.^2-(abs(512-i)).^2);
        fitRangeIndx =512-round(lchord)+1:1:512+round(lchord);
        x = colIndx(fitRangeIndx);
        y = ASItemp(i,fitRangeIndx);
        if length(fitRangeIndx)>10
            pRow(i,:)=polyfit(x,y,nPoly);
        end
        muRow(i,:)=mean(y);
        stdRow(i,:) = std(y);
        a(i).fitRangeIndx = fitRangeIndx;
        background(i,fitRangeIndx) = polyval(pRow(i,:),colIndx(fitRangeIndx));
        if ifit==1
            brightStarIndx = ASItemp(i,:) > (background(i,:)+1*stdRow(i,:));
            if length(find(brightStarIndx)>0)
                ASItemp(i,brightStarIndx) = nan;
                ASItemp(i,fitRangeIndx) = interp_nans(ASItemp(i,fitRangeIndx)')';
            end
        else
            ASINew(i,fitRangeIndx) = ASI(i,fitRangeIndx)-background(i,fitRangeIndx);
            ASINoise(i,fitRangeIndx) = ASItemp(i,fitRangeIndx)-background(i,fitRangeIndx);
        end
        
    end
end

end

function f = fpoly3(p,x)
    f = p(1).*x.^3 + p(2).*x.^2 + p(3).*x + p(4);
end

function [starImage] = faint_star_extracter(ASI, sigma_n)
    
    mSz = 7;
    MSz = 9;
    dSz = (MSz-mSz)/2;
    if mod(dSz,1)
        error('mSz-MSz has to be even');
    end
    
    starImage = zeros(size(ASI));
    
    multiWaitbar('Faint Star Extraction...',0);
    id = 1./size(ASI,1);
     % For faint stars
    for i=1:1:size(ASI,1) - MSz
        for j = 1:1:size(ASI,2) - MSz
        M = ASI(i:i+MSz-1,j:j+MSz-1);
        in = M(1+dSz:dSz+mSz,1+dSz:dSz+mSz);
        out = M;
        out(1+dSz:dSz+mSz,1+dSz:dSz+mSz)=0;
        sig_n = sigma_n(i+dSz+ceil(mSz/2),j+dSz+ceil(mSz/2));
        mu_in = nanmean(in(:));
        mu_out = nanmean(out(:));
        % Condition
        if mu_in > 2*sig_n && mu_out < 2*sig_n 
            starImage(i:i+MSz-1,j:j+MSz-1) = M;
        end
        
        end
    multiWaitbar('Faint Star Extraction...','Increment',id);    
    end

end

function [starImage, imstar] = bright_star_extracter(ASI, sigma_n)
    
    mSz = 11;
    MSz = 13;
    dSz = (MSz-mSz)/2;
    if mod(dSz,1)
        error('mSz-MSz has to be even');
    end
    
    starImage = zeros(size(ASI));
    
    multiWaitbar('Bright Star Extraction...',0);
    id = 1./size(ASI,1);
    imstar = [];
%     countStar = 0;
     % For bright stars
    for i=1:1:size(ASI,1) - MSz
        for j = 1:1:size(ASI,2) - MSz
        M = ASI(i:i+MSz-1,j:j+MSz-1);
        in = M(1+dSz:dSz+mSz,1+dSz:dSz+mSz);
        out = M;
        out(1+dSz:dSz+mSz,1+dSz:dSz+mSz)=0;
        sig_n = sigma_n(i+dSz+ceil(mSz/2),j+dSz+ceil(mSz/2));
        mu_in = nanmean(in(:));
        mu_out = nanmean(out(:));
        % Condition
        if mu_in > 3*sig_n && mu_out < 3*mu_in 
            starImage(i:i+MSz-1,j:j+MSz-1) = M;
%             countStar = countStar + 1;
%             imstar.ID(countStar) = countStar;
%             imstar.image(countStar).M = M;
%             imstar.brightness(countStar) = sum(M(:));
%             imstar.x(countStar) = i;
%             imstar.y(countStar) = j;
        end
        
        end
    multiWaitbar('Bright Star Extraction...','Increment',id);    
    end

end

function imstar = extract_stars(image)
%% Image should be post all processing
binaryImage = zeros(size(image));
binaryImage(image>0)=1;
[labelImage, numSpots] = bwlabel(binaryImage);
props = regionprops(labelImage,image,'Centroid','Area','MeanIntensity');
for i =1:1:numSpots
    imstar.ID = i;
    imstar.location(i,:) = props(i).Centroid;
    imstar.brightness(i) = props(i).Area.*props(i).MeanIntensity;
end
imstar.size = size(image);
end

function newstar = astrometry_file(imstar, el, fileStr)
    if nargin<3
        fileStr = [];
    end
    FOV  = (90-el)*2;
    %% Assuming fish eye
    imLength = imstar.size(1);
    LengthperElevation = imLength/180;
    p.min = round(imLength/2) - round(LengthperElevation*FOV/2);
    p.max = round(imLength/2) + round(LengthperElevation*FOV/2);
    
    selectedStarIndx =(imstar.location(:,1)>p.min &...
        imstar.location(:,1)<p.max & ...
        imstar.location(:,2)>p.min & ...
        imstar.location(:,2)<p.max);
    
    newstar.location = imstar.location(selectedStarIndx,:);
    newstar.brightness = imstar.brightness(selectedStarIndx);
    
    [newstar, I] = sort_star(newstar);

    
    if ~isempty(fileStr)
    dlmwrite(fileStr,round(newstar.location));
    end
    
end

function [newstar, I] = sort_star(newstar)
    [newstar.brightness, I] = sort(newstar.brightness,'descend');
    newstar.brightness = newstar.brightness';
    newstar.location = newstar.location(I,:);
end

function ASI = remove_noise_spikes(ASI, sigma_n)
    
    for i=3:1:size(ASI,1) - 2
        for j = 3:1:size(ASI,2) - 2
            A = ASI(i,j);   
            B1 = ASI(i,j-1);    B2 = ASI(i,j-2);
            C1 = ASI(i-1,j);    C2 = ASI(i-2,j);
            D1 = ASI(i,j+1);    D2 = ASI(i,j+2);
            E1 = ASI(i+1,j);    E2 = ASI(i+2,j);
            s = sigma_n(i,j);
            singleSpike = (A > 2*s)*(B1 < 2*s)*(C1 < 2*s)*(D1 < 2*s)*(E1 < 2*s);
            doubleSpike(1) = (A > 3*s)*((B1>2*s)*(B2<2*s))*...
                not((C1>2*s))*not((D1>2*s))*...
                not((E1>2*s));
            doubleSpike(2) = (A > 3*s)*not((B1>2*s))*...
                ((C1>2*s)*(C2<2*s))*not((D1>2*s))*...
                not((E1>2*s));
            doubleSpike(3) = (A > 3*s)*not((B1>2*s))*...
                not((C1>2*s))*((D1>2*s)*(D2<2*s))*...
                not((E1>2*s));
            doubleSpike(4) = (A > 3*s)*not((B1>2*s))*...
                not((C1>2*s))*not((D1>2*s))*...
                ((E1>2*s)*(E2<2*s));
            ASI(i,j) = not(singleSpike).*A;
            ASI(i,j) = not(sum(doubleSpike)).*A;
            ASI(i,j-1) = not(doubleSpike(1)).*B1;
            ASI(i-1,j) = not(doubleSpike(2)).*C1;
            ASI(i,j+1) = not(doubleSpike(3)).*D1;
            ASI(i+1,j) = not(doubleSpike(4)).*E1;
        end
    end


end











