%% Calibrate Poker Flat ASC
clear all;

stars=get_star_catalogue;
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-full_v2.h5';

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
[polaris.azObs,polaris.elObs] = get_star_az_el...
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
dasc.az = (h5read(fileStr,'/DASC/az'))';
dasc.el = (h5read(fileStr,'/DASC/el'))';
%% Step 2: Identify Hot Pixles
% Load fits File
asiPath = 'C:\Users\nithin\Documents\GitHub\LargeFiles\DASC\20080326\';
asi1File = 'PKR_DASC_0000_20080326_103958.000.FITS';
asi2File = 'PKR_DASC_0000_20080326_104018.000.FITS';
asi3File = 'PKR_DASC_0000_20080326_133018.000.FITS';

ASI1 = fitsread([asiPath,asi1File]);
ASI2 = fitsread([asiPath,asi2File]);
ASI3 = fitsread([asiPath,asi3File]);


[hotPixels] = identify_hot_pixels(ASI1, ASI3, 10);

%% Step 3: Remove background
[ASINew,backgroundRow, ASINoise, pRow, muRow, a] = remove_background(ASI1,5);
[ASINew1,background1, ASINoise1] = remove_background(ASINew');
ASINew1 = ASINew1';
ASINoise1 = ASINoise1';
background1 = background1';

figure; 
imagesc(ASINew);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed (Row fit)');
figure; 
imagesc(ASINoise);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 10]); title('ASI noise (Row fit)');
figure;
imagesc(background);
colorbar;
colormap(get_colormap('k','w'));
caxis([200 800]); title('Background with stars removed (Row fits)');
%%
figure; 
imagesc(ASINew1);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed (Col fit)');
%%
figure;
imagesc(ASINoise1);
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
totalNoise = ASINoise1+ASINoise;
totalNoise(totalNoise==0)=nan;
sigma_n = nlfilter(totalNoise,[9 9],fnoise);
%%
figure; 
imagesc(sigma_n);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 10]); title('\sigma_n Noise at each pixel');

%% Step 5: Noise spike removal
ASINew2 = remove_noise_spikes(ASINew1, sigma_n);
figure; 
imagesc(ASINew2);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('ASI with bkgnd removed + noise spikes removed (Col fit)');

%% Step 6: Star extraction 
faintStarImage = faint_star_extracter(ASINew1, sigma_n);
% Rough removal of noise below 20 pixels
%
faintStarImage(faintStarImage <= 20) = 0;
faintStarImage(1:125,:) = 0;
faintStarImage(1000:1024,:) = 0;

%%
brightStarImage = bright_star_extracter(ASINew1, sigma_n);
% Rough removal of noise below 20 pixels
%
brightStarImage(brightStarImage <= 20) = 0;
brightStarImage(1:125,:) = 0;
brightStarImage(1000:1024,:) = 0;

%%
figure; 
imagesc(faintStarImage);
colorbar;
colormap(get_colormap('k','w'));
caxis([0 100]); title('Faint stars (extraction algorithm)');
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
plot(x,ASI1(i,x));
% ylim([0,1000]);
%% Plotting
figure;
colormap(viridis);
% Plotting all-sky camera image
indx = dasc.el>0;
plot_DASC_aer(dasc.ASI(indx), dasc.az(indx), dasc.el(indx), 512);
hold on; 
plot_grid_aer([0,90],[22.5],'m');
% colorbar;
xlim([-120,+120]);
ylim([-120,+120]);
caxis([300 450]);
hold on;

% Stars
starfilter=stars.vmag<3;
plot_aer_stars(stars.az(starfilter), stars.el(starfilter), stars.relIntensity(starfilter)*50,'r');
hold on;
plot_aer_stars(polaris.azStellaris,polaris.elStellaris,10,'w');
hold on;
plot_aer_star_label(polaris.az, polaris.el, 'Polaris');
hold on;
plot_aer_stars(polaris.azObs, polaris.elObs,10,'c');

title(timeStr);

function plot_aer_stars(az,el,relIntensity,colorStr)
    x = (90-el).*sind(360-az);
    y = (90-el).*cosd(az);
    
    scatter(x,y,relIntensity,colorStr);
    
end

function plot_aer_star_label(az,el,textString,colorStr)
    x = (90-el).*sind(360-az);
    y = (90-el).*cosd(az);
    
    if nargin<4
        colorStr = 'r';
    end
    
    text(x, y, textString, 'Color', colorStr);
    
end

function [totalIntensity,timeArr]=estimate_darkest_frame(h5FileStr)

asi = permute(h5read(h5FileStr,'/DASC/ASI'),[3 2 1]);
timeArr = unixtime2matlab((h5read(h5FileStr,'/DASC/time'))');
totalIntensity = sum(sum(asi,3),2);

end

function [hotPixels] = identify_hot_pixels(ASI1, ASI2, threshold)
    
    if nargin<3
        threshold = 2;
    end
    
    temp = ASI1 - ASI2;
    temp2 = temp;
    hotPixels=zeros(size(temp2));
    hotPixels(abs(temp2)<=threshold) = 1;
    
end

function [ASINew, background, ASINoise, pRow, muRow, a] = remove_background(ASI, nPoly)

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

function starImage = faint_star_extracter(ASI, sigma_n)
    
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
        if mu_in > 3*sig_n && mu_out < 2*sig_n 
            starImage(i:i+MSz-1,j:j+MSz-1) = M;
        end
        
        end
    multiWaitbar('Faint Star Extraction...','Increment',id);    
    end

end

function starImage = bright_star_extracter(ASI, sigma_n)
    
    mSz = 11;
    MSz = 13;
    dSz = (MSz-mSz)/2;
    if mod(dSz,1)
        error('mSz-MSz has to be even');
    end
    
    starImage = zeros(size(ASI));
    
    multiWaitbar('Bright Star Extraction...',0);
    id = 1./size(ASI,1);
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
        if mu_in > 10*sig_n && mu_out < 0.1*mu_in 
            starImage(i:i+MSz-1,j:j+MSz-1) = M;
        end
        
        end
    multiWaitbar('Bright Star Extraction...','Increment',id);    
    end

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











