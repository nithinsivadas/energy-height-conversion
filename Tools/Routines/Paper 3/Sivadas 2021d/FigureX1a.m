%% ACF of the solar wind Electric field 
% In this routine 
% - we find the ACF of solar wind electric field from wind measurements
% - we device a way to construct a stochastic process
% with an underlying lognormal distribution, with time series that 
% has an ACF equivalent to the solar wind electric field. 

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

%%
Vmag = wind1.velocity; 
Bmag = (wind1.BxGSM.^2+wind1.ByGSM.^2 + wind1.BzGSM.^2).^0.5; 

% filter.geo = ((geotail.noseYGSE-geotail.YGSE).^2 + (geotail.noseZGSE-geotail.ZGSE).^2).^0.5 < 20;

% Data from spacecraft
% XS = geotail.E_kl(filter.geo)*10^-3;
% XS = wind1.E_kl(~isnan(wind1.E_kl))*10^-3;

% XS1 = wind1.velocity(~isnan(wind1.velocity));
% ff = ~isnan(wind1.BxGSM) & ~isnan(wind1.ByGSM) & ~isnan(wind1.BzGSM);

% XS2 = Bmag(ff);
% XS = XS1 - nanmean(XS1);
% XSS = XS2 - nanmean(XS2);

% dt = 1; %minute
% lag = fftshift(-length(XS):1:length(XS)-1)'.*dt;
% [RArray ,lag1]= xcorr(XS,'coeff');
% [RArray2 ,lag2]= xcorr(XSS,'coeff');
% RArray(2:end+1,:) = RArray(1:end,:);
% RArray(1,:)=RArray(2,:);
% RArray = ifftshift(RArray);
% plot(lag1/60,RArray);
% xlim([0,80]);
% hold on;
% plot(lag2/60,RArray2);


%%
[Mv, nSample,  nEnsemble]  = split_series(Vmag,2^13);
[RArrayv,lagv] = find_correlation(Mv - nanmean(Mv,1),nSample,nEnsemble,1);
%%
[Mb, nSample,  nEnsemble]  = split_series(Bmag,2^13);
[RArrayb,lagb] = find_correlation(Mb-nanmean(Bmag,1),nSample,nEnsemble,1);
%%
[Me, nSample,  nEnsemble]  = split_series(wind1.E_kl.*10^-3,2^13);
[RArraye,lage] = find_correlation(Me-nanmean(Me,1),nSample,nEnsemble,1);
%%
figure;
plot(fftshift(lagv/60), fftshift(mean(RArrayv./RArrayv(:,1))));
hold on;
plot(fftshift(lagb/60), fftshift(mean(RArrayb./RArrayb(:,1))));

hold on;
plot(fftshift(lage/60), fftshift(mean(RArraye./RArraye(:,1))));
xlim([0,80]);
legend('V_{sw}','B_{mag}','E_{M}');
%%
acf_E = mean(RArraye(:,1:2^13)./RArraye(:,1));
lag_E = lage(1:2^13);
f = fit(lag_E,acf_E','spline'); % Rxx(t) = 0.4968*exp(-0.2874*t) + 0.2911*exp(-0.02941*t) [hrs]
% Rxx(t) = 0.4872*exp(-0.005368*t) + 0.3142*exp(-0.0005358*t) [mins]
figure; 
plot(lag_E/60,log(acf_E));
hold on;
plot(lag_E/60,log(f(lag_E)));
set(gca,'XScale','linear');

%% Doing the cholensky decomposition to get a random variable that correrlates
dt = 1; 
nSamples = 2^13;
nEnsemble=1000;
pdfx1 = makedist('Lognormal',-0.286163,1.09028);
% pdfx1 = makedist('Lognormal',-0.286163,1.09028);
% pdfx1 = makedist('Normal',0,1.09028);
X1 = random(pdfx1,nSamples,nEnsemble);


% X = random(pdfx,nSamples,1);
lag = fftshift(-nSamples:1:nSamples-1)'.*dt;
% Rm = 0.4968*exp(-0.2874*abs(lag)/60) + 0.2911*exp(-0.02941*abs(lag)/60);
% Rm = 0.4968*exp(-0.2874*abs(lag)*60) + 0.2911*exp(-0.02941*abs(lag)*60);
% Rm = 0.4872*exp(-0.005368*abs(lag)) + 0.3142*exp(-0.0005358*abs(lag));
Rm = f(abs(lag));
% Rm = 1*exp(-0.005368*abs(lag)) ;
RmMatrix = toeplitz(Rm(find(lag==0):find(lag==(nSamples-1)*dt)));


Y = MvLogNRand(repmat(-0.286163,nSamples,1),repmat(1.09028,nSamples,1),nEnsemble,RmMatrix);

% B = chol(RmMatrix,'lower');
% Y = mean(X1) + (X1-mean(X1))'*B; 
%%
figure;
XBin = logspace(-2,3,100);
histogram(X1(:),XBin,'Normalization','pdf');
hold on;
histogram(Y(:),XBin,'Normalization','pdf');
set(gca,'XScale','log');
%%
[Rt,lagt]=find_correlation(Y-mean(Y,1),nSamples,nEnsemble,1);
%%
figure; 
plot(fftshift(lagt/60),fftshift(mean(Rt./Rt(:,1))));
xlim([-100,100]);
hold on;
plot(fftshift(lag/60),fftshift(Rm));


function [M, nSample, nEnsemble]= split_series(series,sampleSize,MaxNoOfMissingValues)
    if nargin<3
        MaxNoOfMissingValues=1000;
    end
        
    series = padarray(series,sampleSize-mod(length(series),sampleSize),'post');
    L = length(series); 
    M0 = reshape(series,sampleSize,[])';
    indx=sum(isnan(M0),2)>MaxNoOfMissingValues;
    l = 1:1:size(M0,1);
    k=1;
    for i=l(~indx)
        M(k,:) = M0(i,:);
        k=k+1;
    end
    nSample = size(M,2);
    nEnsemble = size(M,1);
    
    M = interp_nans(M')';
    
%     M = M-mean(M,2);
end

function [RArray,lag]=find_correlation(M,nSample,nEnsemble,sampleTime)
    dt = sampleTime;
    lag = fftshift(-nSample:1:nSample-1)'.*dt;
    xCell = mat2cell(M, ones(1,nEnsemble),nSample);
    [RArray] = cellfun(@(x) xcorr(x,'unbiased'),xCell,'UniformOutput',false);
    RArray = cell2mat(RArray);
    RArray(:,2:end+1) = RArray(:,1:end);
    RArray(:,1)=RArray(:,2);
    RArray = ifftshift(RArray);
end

function y = MvLogNRand( Mu , Sigma , Simulations , CorrMat )
%MVLOGNRAND MultiVariant Lognormal random numbers with correlation
%
%   Mu: The Lognormal parameter Mu  (can be column or row vector)
%
%   Sigma: The Lognormal parameter Sigma (can be column or row vector)
%
%   Simulations:  The Number of simulations to run (scalar)
%
%   CorrMat:  OPTIONAL A square matrix with the number of rows and columns
%   equal to the number of elements in Mu/Sigma.  Each element on the
%   diagonal is equal to one, with the off diagonal cells equal to the
%   correlation of the marginal Lognormal distributions. If not specified,
%   then assume zero correlation.
%
%   To check the simulation run corrcoef(Y) and that should be the same as
%   your CorrMat.
%
%   REQUIRES THE STATISTICS TOOLBOX
%
%   Example:
%   Mu    = [ 11 12 13 ];
%   Sigma = [ .1 .3 .5 ];
%   Simulations = 1e6;
%   CorrMat = [1 .2 .4 ; .2 1 .5 ; .4  .5 1];
%   y = MvLogNRand( Mu , Sigma , Simulations , CorrMat );
%
%   corrcoef(y)
%   ans =
%            1      0.19927      0.40156
%      0.19927            1      0.50008
%      0.40156      0.50008            1
%
%   CorrMat =
%               1          0.2          0.4
%             0.2            1          0.5
%             0.4          0.5            1
%
%   For more information see: Aggregration of Correlated Risk Portfolios:
%   Models and Algorithms; Shaun S. Wang, Phd.  Casualty Actuarial Society
%   Proceedings Volume LXXXV www.casact.org
%
%   Author: Stephen Lienhard
% Error checking
if nargin < 3
    error('Must have at least 3 input arguements')
end
if numel(Simulations) ~= 1 || Simulations < 0
    error('The number of simulations must be greater then zero and a scalar')
end
if nargin == 3
    CorrMat = eye(numel(Mu));
elseif size(CorrMat,1) ~= size(CorrMat,2)
    error('The correlation matrix must have the same number of rows as columns')
end
if numel(Mu) ~= numel(Sigma)
    error('Mu and Sigma must have the same number of elements')
end
% Force column vectors
Mu     = Mu(:);
Sigma  = Sigma(:);
% Calculate the covariance structure
sigma_down = repmat( Sigma' , numel(Sigma), 1            );
sigma_acrs = repmat( Sigma  , 1           , numel(Sigma) );
covv = log( CorrMat .* sqrt(exp(sigma_down.^2)-1) .* ...
                       sqrt(exp(sigma_acrs.^2)-1) + 1 );
% The Simulation
y = exp( mvnrnd( Mu , covv , Simulations ));
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
