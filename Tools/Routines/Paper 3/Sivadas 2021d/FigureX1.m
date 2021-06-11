%% Figure X1
%% Plot all errors, distributions and its effect
% A toy-model,based on data. 
% Date: 8th June 2021

%%
omniFile = 'G:\My Drive\Research\Projects\Data\omni.h5';
omni = extract_omni_data(omniFile);


%% Extract error and X pdf from data
dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';
geotail1 = extract_set_of_data(dataFolder,'geotail');
wind1 = extract_set_of_data(dataFolder,'wind');
ace1 = extract_set_of_data(dataFolder,'ace');

% ACE, WIND, GEOTAIL, time conjunctions
[~,~,wi]=intersect(geotail1.datetime,wind1.datetime);
[~,ai,wii]=intersect(ace1.datetime,wind1.datetime(wi));
[~,gi,~]=intersect(geotail1.datetime,wind1.datetime(wi(wii)));
wi = wi(wii);

% ACE, WIND, GEOTAIL conjunction tables
ace = ace1(ai,:);
wind = wind1(wi,:);
geotail = geotail1(gi,:);

%% Calculate error and X pdf from data

% Filter the data to be within 20 RE of the nose
filter.geo = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 20;

% Data from spacecraft
XS = geotail.E_kl(filter.geo)*10^-3;
XmSa = wind.E_kl(filter.geo)*10^-3; % Measured X 
XmSb = ace.E_kl(filter.geo)*10^-3;
% Estimating WIND/ACE measurements, if there is a delay of 8 min from true
% value
sigmaT = 8;
time = geotail.datetime(filter.geo);
time_Sa_N = time + minutes(random('Weibull',sigmaT,2,length(time),1));
XmSaT_Sa_N = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(time_Sa_N));
% Estimating Geotail-Nose delay, if there is a delay 1 min from true value
sigmaT = 1;
time_S_N = time + minutes(random('Weibull',sigmaT,2,length(time),1));
XST_S_N = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(time_S_N));
% Estimating Nose-Ionosphere delay, if there is a delay 20 min from true value
sigmaT = 20;
time_N_I = time + minutes(random('Weibull',sigmaT,2,length(time),1));
XST_N_I = interp1(datenum(geotail1.datetime),geotail1.E_kl.*10^-3,datenum(time_N_I));


% Error
eS = XmSa - XS; 
eS2 = XmSa - XmSb;
% Errors due to time delay

eST = XmSaT_Sa_N - XS; % Artifically constructed XmSa, through 8 min propagation delay error
eS_S_N = XST_S_N - XS; % Additional error on true value, cause of the propagation delay with the nose
eS_N_I = XST_N_I - XS; % Additional error on true value, cause of propagation delay from nose to ionosphere

eTot = eST+eS_N_I; %Total error, for now just error with time delay 8 min
eNormSTot = eTot./XS;

eNormS =  eS./XS;
eNormS2 =  eS2./XS;

% Identifying the kernel distributions
% pdfe = fitdist(eNormS(~isnan(eNormS) & XS>0.1 & XS<30),'kernel'); % Normalized error distribution
pdfeTotal = fitdist(eNormSTot(~isnan(eNormSTot) & XS>2 & XS<30),'kernel'); % Normalized error distribution for just 8 min delay
% pdfx = fitdist(XS(~isnan(XS) & XS>0 & XS<30),'kernel');
% pdfx = makedist('Lognormal',-0.286163,0.789428);
% Calculating the Toy model
%
Nsamples = 1051200;
Nensemble = 1;

X = random(pdfx,Nsamples,Nensemble); % True values of variable of interest



time1 = 1:1:length(X);
sigmaT = 8;
time_a_N = time1 + (random('Weibull',sigmaT,2,length(time1),1))';
X_a_N = interp1(time1,X,time_a_N)';
e_a_N = X_a_N - X;
eNorm_a_N = e_a_N./X; 
pdfe = fitdist(eNorm_a_N,'kernel');

U = random(pdfe,Nsamples,Nensemble); % Normalized error (almost homoscedastic)
UT = random(pdfeTotal,Nsamples,Nensemble);

W = X + X.*U; % Errorneous measurements made by WIND-like spacecraft
Wa = X + X.*random(pdfe,Nsamples,Nensemble); % Another spacecraft measurement (like ACE)

W2 = X + X.*UT; % Total error, but right now just from 8 min uncertainty

% Error
e = W - X; 
e2 = W - Wa;
eT = W2 - X;

% Normalized error
eNorm = e./X; % Type 1 - difference between errenous measurements and true-value
eNorm2 = e2./X; % Type 2 - difference between two errenous measurements
eNormT = eT./X; 

% Calculate expectations and PDF
[PegX,XPegX,YPegX] = conditional_pdf(e,X,-5:0.1:5,0:1:20);
[PeTgX,XPeTgX,YPeTgX] = conditional_pdf(eT,X,-5:0.1:5,0:1:20); %Total/time

[PeNgX,XPeNgX,YPeNgX] = conditional_pdf(eNorm,X,-5:0.1:5,0:1:20);
[PeNTgX,XPeNTgX,YPeNTgX] = conditional_pdf(eNormT,X,-5:0.1:5,0:1:20); %Total/time

XBins = 0:1:30;

EXgW = create_curve(W,X,XBins);
EXgW2 = create_curve(W2,X,XBins);
EXgX = create_curve(X,X,XBins);

%
YBins = -100:-100:-3000;
Y = -(95.*X + random('Normal',0,1,Nsamples,Nensemble));
Yw = Y + random('Normal',0,2,Nsamples,Nensemble);
YS = omni.Fsml(datenum(time));
EYgW = create_curve(W,Yw,XBins);
EYgW2 = create_curve(W2,Yw,XBins);
EYgX = create_curve(X,Yw,XBins);
EYSgXm = create_curve(XmSa,YS,XBins);

% Figure
set(0,'defaulttextInterpreter','latex');

f = figure;
resize_figure(f);
p = panel();
p.marginright=20;
p.pack(4,2);

p.select('all');


XBins1 = logspace(-2,2,100);
eBins = -5:0.1:5;
p(1,1).select();
histogram(X,XBins1,'Normalization','pdf'); 
hold on;
histogram(W,XBins1,'Normalization','pdf'); 
set(gca,'XScale','log','YScale','log','XTick',[0.1,1,10,30]);
legend('True','Measured');
xlabel('X');
ylabel('pdf');

% ylim([10^-5, 1]);
p(1,2).select();
histogram(U,eBins,'Normalization','pdf');
hold on;
histogram(UT,eBins,'Normalization','pdf');
xlabel('$(X^*-X)/X$','Interpreter','latex');
ylabel('pdf');
legend('Spatial','Temporal (8 min)');

p(2,1).select();
colormap(inferno);
plot_2D_error(XPegX, YPegX, PegX);
caxis([0,1]);
xlabel('$P(X^*-X|X)$','Interpreter','latex');
ylabel('$X$');
title('Spatial Uncertainty');

p(2,2).select();
colormap(inferno);
plot_2D_error(XPeTgX, YPeTgX, PeTgX);
caxis([0,1]);
xlabel('$P(X^*-X|X)$','Interpreter','latex');
ylabel('$X$');
title('Temporal Uncertainty (8 min propagation delay)');

p(3,1).select();
colormap(inferno);
plot_2D_error(XPeNgX, YPeNgX, PeNgX);
caxis([0,1]);
xlabel('$P((X^*-X)/X|X)$','Interpreter','latex');
ylabel('$X$');
title('Spatial Uncertainty');
hold on;
plot(repmat(0,100),linspace(0,20,100),'k');

p(3,2).select();
colormap(inferno);
plot_2D_error(XPeNTgX, YPeNTgX, PeNTgX);
caxis([0,1]);
xlabel('$P((X^*-X)/X|X)$','Interpreter','latex');
ylabel('$X$');
title('Temporal Uncertainty (8 min propagation delay)');
hold on;
plot(repmat(0,100),linspace(0,20,100),'k');

p(4,1).select();
p1 = plot_curve(EXgX,'k');
hold on;
p2 = plot_curve(EXgW,'b');
hold on;
p3 = plot_curve(EXgW2,'r');

p(4,2).select();
p1 = plot_curve(EYgX,'k');
hold on;
p2 = plot_curve(EYgW,'b');
hold on;
p3 = plot_curve(EYgW2,'r');
hold on;
p4 = plot_curve(EYSgXm,'g');


figure; 
curve1 = create_curve(XmSa,eTot,0:0.5:30);
plot(curve1.XBins,100.*curve1.stdYgX./curve1.XBins);
hold on;
curve2 = create_curve(W2,eT,0:0.5:30);
plot(curve2.XBins,100.*curve2.stdYgX./curve2.XBins);
legend('Satellite Data','Toy model');
xlabel('$X^*$');
ylabel('$std(X^*)/<X^*>$');
%% Functions

function RmArray = autocorrelation_model(N,A,fWidth,omegaDoppler,lag)
       RmArray =N*double(lag==0) + A.*exp(-abs(lag).*fWidth).*exp(-1i*omegaDoppler.*lag);
end
   

function plot_2D_error(Y,X,P)
    p = pcolor(Y,X,P);
    set(p,'EdgeColor','none');
    colorbar_thin();
end

function curve = create_curve(X, Y, Ei)

if nargin<3
    Ei = 100;
end
Y = Y(:);
Y(Y==999999)=nan;

X1 = X(~isnan(X) & ~isnan(Y));
Y1 = Y(~isnan(X) & ~isnan(Y));
X = X1;
Y = Y1;

[xindx, E] = discretize(X(:),Ei);

for i = 1:max(xindx)
    curve.YgX(i) = nanmean(Y(xindx==i));
    curve.stdYgX(i) = nanstd(Y(xindx==i));
    curve.NSamples(i) = sum(xindx==i & ~isnan(Y));
    curve.SEM(i) = nanstd(Y(xindx==i))./sqrt(curve.NSamples(i));
    curve.ts(i,:) = tinv([0.025 0.975],curve.NSamples(i)-1);
    curve.CI(i,:) = curve.YgX(i) + curve.ts(i,:)*curve.SEM(i);
    curve.XBins(i) = 0.5*(E(i)+E(i+1));
end

curve.E = E; 

end

function [fXgY,XX,YY] = conditional_pdf(X, Y, gridx, gridy)

X1 = X(~isnan(X) & ~isnan(Y));
Y1 = Y(~isnan(X) & ~isnan(Y));
X = X1;
Y = Y1;

sz = 2^10;
Y = Y(:);
Y(Y==999999)=nan;

[bandwidth,fXY,XX,YY]=kde2d([X,Y],sz,[min(gridx),min(gridy)],[max(gridx),max(gridy)]);
[fY, YY1] = ksdensity(Y,YY(:,1));

fXgY = fXY./repmat(fY,1,sz);

% [XX, YY] = meshgrid(gridx,gridy);
% 
% [fxy, yii] = ksdensity([X,Y],[XX(:),YY(:)]);
% Fxy = scatteredInterpolant(yii(:,1),yii(:,2),fxy);
% 
% [fy, yi] = ksdensity(Y,gridy);
% Fy = griddedInterpolant(yi,fy);
% 
% cpdf = Fxy(XX,YY)./(Fy(gridy)')./(sum((Fxy(XX,YY)./(Fy(gridy)')),2));

end

function p = plot_curve(curve, color)

CI1 = interp_nans(curve.CI);
p=plot(curve.XBins, curve.YgX,color);
hold on;
plot_ci(curve.XBins,CI1,color,0.2);
end

function plot_ci(x,ci,color,alpha)
    hold on;
    X2 = [x, fliplr(x)];
    inBetween = [ci(:,1)', fliplr(ci(:,2)')];
    fill(X2,inBetween,color,'LineStyle','none','FaceAlpha',alpha);
end

function T = extract_set_of_data(loadFolder,satellite)
    fileStr = get_files_in_folder(loadFolder,[satellite,'_min_b*.txt']);
    for i=1:1:length(fileStr)
        T1 = extract_data([loadFolder,fileStr{i}]);
        if i==1
            T=T1;
        else
            T = [T;T1];
        end
    end
end

function [bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
% fast and accurate state-of-the-art
% bivariate kernel density estimator
% with diagonal bandwidth matrix.
% The kernel is assumed to be Gaussian.
% The two bandwidth parameters are
% chosen optimally without ever
% using/assuming a parametric model for the data or any "rules of thumb".
% Unlike many other procedures, this one
% is immune to accuracy failures in the estimation of
% multimodal densities with widely separated modes (see examples).
% INPUTS: data - an N by 2 array with continuous data
%            n - size of the n by n grid over which the density is computed
%                n has to be a power of 2, otherwise n=2^ceil(log2(n));
%                the default value is 2^8;
% MIN_XY,MAX_XY- limits of the bounding box over which the density is computed;
%                the format is:
%                MIN_XY=[lower_Xlim,lower_Ylim]
%                MAX_XY=[upper_Xlim,upper_Ylim].
%                The dafault limits are computed as:
%                MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
%                MAX_XY=MAX+Range/4; MIN_XY=MIN-Range/4;
% OUTPUT: bandwidth - a row vector with the two optimal
%                     bandwidths for a bivaroate Gaussian kernel;
%                     the format is:
%                     bandwidth=[bandwidth_X, bandwidth_Y];
%          density  - an n by n matrix containing the density values over the n by n grid;
%                     density is not computed unless the function is asked for such an output;
%              X,Y  - the meshgrid over which the variable "density" has been computed;
%                     the intended usage is as follows:
%                     surf(X,Y,density)
% Example (simple Gaussian mixture)
% clear all
%   % generate a Gaussian mixture with distant modes
%   data=[randn(500,2);
%       randn(500,1)+3.5, randn(500,1);];
%   % call the routine
%     [bandwidth,density,X,Y]=kde2d(data);
%   % plot the data and the density estimate
%     contour3(X,Y,density,50), hold on
%     plot(data(:,1),data(:,2),'r.','MarkerSize',5)
%
% Example (Gaussian mixture with distant modes):
%
% clear all
%  % generate a Gaussian mixture with distant modes
%  data=[randn(100,1), randn(100,1)/4;
%      randn(100,1)+18, randn(100,1);
%      randn(100,1)+15, randn(100,1)/2-18;];
%  % call the routine
%    [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none'), view([0,60])
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
%
% Example (Sinusoidal density):
%
% clear all
%   X=rand(1000,1); Y=sin(X*10*pi)+randn(size(X))/3; data=[X,Y];
%  % apply routine
%  [bandwidth,density,X,Y]=kde2d(data);
%  % plot the data and the density estimate
%  surf(X,Y,density,'LineStyle','none'), view([0,70])
%  colormap hot, hold on, alpha(.8)
%  set(gca, 'color', 'blue');
%  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
%
%  Reference:
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
global N A2 I
if nargin<2
    n=2^8;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;
N=size(data,1);
if nargin<3
    MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
    MAX_XY=MAX+Range/2; MIN_XY=MIN-Range/2;
end
scaling=MAX_XY-MIN_XY;
if N<=size(data,2)
    error('data has to be an N by 2 array where each row represents a two dimensional observation')
end
transformed_data=(data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
%bin the data uniformly using regular grid;
initial_data=ndhist(transformed_data,n);
% discrete cosine transform of initial data
a= dct2d(initial_data);
% now compute the optimal bandwidth^2
  I=(0:n-1).^2; A2=a.^2;
 t_star=root(@(t)(t-evolve(t)),N);
p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
% smooth the discrete cosine transform of initial data using t_star
a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
% now apply the inverse discrete cosine transform
if nargout>1
    density=idct2d(a_t)*(numel(a_t)/prod(scaling));
	density(density<0)=eps; % remove any negative density values
    [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
end
bandwidth=sqrt([t_x,t_y]).*scaling; 
end
%#######################################
function  [out,time]=evolve(t)
global N
Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
time=(2*pi*N*Sum_func)^(-1/3);
out=(t-time)/time;
end
%#######################################
function out=func(s,t)
global N
if sum(s)<=4
    Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
    time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
    out=psi(s,time);
else
    out=psi(s,t);
end
end
%#######################################
function out=psi(s,Time)
global I A2
% s is a vector
w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
wx=w.*(I.^s(1));
wy=w.*(I.^s(2));
out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
end
%#######################################
function out=K(s)
out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
end
%#######################################
function data=dct2d(data)
% computes the 2 dimensional discrete cosine transform of data
% data is an nd cube
[nrows,ncols]= size(data);
if nrows~=ncols
    error('data is not a square array!')
end
% Compute weights to multiply DFT coefficients
w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
weight=w(:,ones(1,ncols));
data=dct1d(dct1d(data)')';
    function transform1d=dct1d(x)
        % Re-order the elements of the columns of x
        x = [ x(1:2:end,:); x(end:-2:2,:) ];
        % Multiply FFT by weights:
        transform1d = real(weight.* fft(x));
    end
end
%#######################################
function data = idct2d(data)
% computes the 2 dimensional inverse discrete cosine transform
[nrows,ncols]=size(data);
% Compute wieghts
w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
weights=w(:,ones(1,ncols));
data=idct1d(idct1d(data)');
    function out=idct1d(x)
        y = real(ifft(weights.*x));
        out = zeros(nrows,ncols);
        out(1:2:nrows,:) = y(1:nrows/2,:);
        out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
    end
end
%#######################################
function binned_data=ndhist(data,M)
% this function computes the histogram
% of an n-dimensional data set;
% 'data' is nrows by n columns
% M is the number of bins used in each dimension
% so that 'binned_data' is a hypercube with
% size length equal to M;
[nrows,ncols]=size(data);
bins=zeros(nrows,ncols);
for i=1:ncols
    [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
    bins(:,i) = min(bins(:,i),M);
end
% Combine the  vectors of 1D bin counts into a grid of nD bin
% counts.
binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t=root(f,N)
% try to find smallest root whenever there is more than one
N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
tol=10^-12+0.01*(N-50)/1000;
flag=0;
while flag==0
    try
        t=fzero(f,[0,tol]);
        flag=1;
    catch
        tol=min(tol*2,.1); % double search interval
    end
    if tol==.1 % if all else fails
        t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
    end
end
end

function omni = extract_omni_data(omniFile)
    
    %Time
    time = unixtime2matlab(h5read(omniFile,'/Time'));
    omni.time = time; 
    
        % Auroral electrojet indices
    SML = h5read(omniFile,'/Indices/SML');
    omni.Fsml = griddedInterpolant(time,SML);
    
    SMU = h5read(omniFile,'/Indices/SMU');
    omni.Fsmu = griddedInterpolant(time,SMU);
    
    AL = h5read(omniFile,'/Indices/AL');
    AL(AL==99999)=nan;
    omni.FAL = griddedInterpolant(time, AL);
    
    AU = h5read(omniFile,'/Indices/AU');
    AU(AU==99999)=nan;
    omni.FAU = griddedInterpolant(time, AU);
    
    symH = h5read(omniFile,'/Indices/SYM_H');
    omni.FsymH = griddedInterpolant(time, symH);
    % Solar wind
    
        % IMF Bz
    BzGSM = h5read(omniFile,'/BField/BzGSM');
    omni.FBz = griddedInterpolant(time,BzGSM);
    
        % IMF By
    ByGSM = h5read(omniFile,'/BField/ByGSM');
    omni.FBy = griddedInterpolant(time,ByGSM);
    
    
        % IMF B_T (Tangential to to GSM_x)
    BzGSM(BzGSM==0) = 0.0001;
    B_T = (ByGSM.^2 + BzGSM.^2).^0.5;
    
        % IMF Clock angle
    theta_kl = wrapTo2Pi(atan2(B_T,BzGSM));
        
        
        % Velocity
    velocity = h5read(omniFile,'/Velocity/V');
    omni.Fv = griddedInterpolant(time,velocity);
    
    
    % Solarwind - Magnetosphere Coupling
    l_0 = 7*(6371*10^3);
    
    E_kl = velocity.*B_T.*(sin(theta_kl/2)).^2; %Km nT/s
    % The “geoeffective” (or “merging”) electric field [Kan and Lee, 1979]
    omni.Fekl = griddedInterpolant(time,E_kl);    
    
end

function T1 = extract_data(loadFile)
 
        format = '%4f %4f %3f %3f %4f %4f %4.1f %6f %6.2f %6.2f %6.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %6.2f %8.2f %8.2f %4f %8.1f %8.1f %8.1f %8.1f %7.2f %9.0f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %7f';
        tempData = load_ascii_files(loadFile, format, 0);
        T.datetime = datetime(tempData{1},1,tempData{2},tempData{3},tempData{4},zeros(size(tempData{4})));
        T.timeshift = tempData{8}; T.timeshift(T.timeshift==999999) = nan;
        T.BxGSM = tempData{13};  T.BxGSM(T.BxGSM > 9999) = nan;
        T.ByGSM = tempData{16}; T.ByGSM(T.ByGSM > 9999) = nan;
        T.BzGSM = tempData{17}; T.BzGSM(T.BzGSM > 9999) = nan;
        T.velocity = tempData{23}; T.velocity(T.velocity > 99999) = nan;
        T.B_T = (T.ByGSM.^2 + T.BzGSM.^2).^0.5;
        T.theta_kl = wrapTo2Pi(atan2(T.B_T,T.BzGSM));
        T.E_kl = T.velocity.*T.B_T.*(sin(T.theta_kl/2)).^2;
        T.XGSE = tempData{29}; T.XGSE(T.XGSE > 9999) = nan;
        T.YGSE = tempData{30}; T.YGSE(T.YGSE > 9999) = nan;
        T.ZGSE = tempData{31}; T.ZGSE(T.ZGSE > 9999) = nan;
        T.noseXGSE = tempData{32}; T.noseXGSE(T.noseXGSE > 9999) = nan;
        T.noseYGSE = tempData{33}; T.noseYGSE(T.noseYGSE > 9999) = nan;
        T.noseZGSE = tempData{34}; T.noseZGSE(T.noseZGSE > 9999) = nan;
        T.distance = sqrt((T.XGSE-T.noseXGSE).^2 + (T.YGSE-T.noseYGSE).^2 + (T.ZGSE-T.noseZGSE).^2);
        T1 = struct2table(T);
end
