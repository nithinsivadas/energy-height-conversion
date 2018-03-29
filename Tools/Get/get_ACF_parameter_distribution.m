%% Monte-cralo estimation of fitted ISR spectra/ACF and the distribution of its parameters
% INTRODUCTORY TEXT
%%
function [noiseDist,signalDist,dopplerShiftDist,fWidthDist]...
    = get_ACF_parameter_distribution(nSamples,nPulseTrains,nExperiments,...
    noisePower,signalPower,dopplerFrequency,fWidth,sampleTime)
%% get_ACF_parameter_distribution
%   Runs a monte-carlo simulation to estimate the statistical distribution of singal,
%   noise, doppler-shifted frequency and spectra width from a modelled stochastic ISR
%   signal with an underlying probability distribution. The probability distribution is 
%   gaussian with zero mean and an autocorrelation function modelled with parameters 
%   specified by noisePower, signalPower, dopplerFrequency, fWidth.
%   $x^2+e^{\pi i}$
%% Initializing 
% Sample Time
dt = sampleTime;

% Lags of the model autocorrelation function
lag = fftshift(-nSamples:1:nSamples-1)'.*dt; %Has 2*nSamples elements for fastest fft

% Autocorrelation model parameters
N = noisePower;
A = signalPower;
omegaDoppler = 2*pi*dopplerFrequency;

%Array of autocorrelation functions for each pulse train and experiment
RArray3D = zeros(nSamples*2,nPulseTrains,nExperiments); 

%Output value of the Non-lineaR least square fitting routine/fitted parameters 
outputNLR = zeros(4,nPulseTrains,nExperiments);

%Initial value of the parameters that we are fitting for
inputNLR0=[0.1,10.1,0.1,1.1]; %Noise, Signal Strength, wDoppler, fWidth

%Fitting options
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',10000,'MaxIterations',2000,...
    'TypicalX',[10,10,10,10],'Display','off');

%% Step2: Generating voltage samples that follow model ACF
% Calculating the model ACF at each lag
Rm = autocorrelation_model(N,A,fWidth,omegaDoppler,lag);
% Calculating the Autocorrelation Matrix
RmMatrix = toeplitz(Rm(find(lag==0):find(lag==(nSamples-1)*dt))).';
% Taking effectively the matrix square-root (Cholesky factorization)
B = chol(RmMatrix,'lower');

for thisExperiment = 1:1:nExperiments
    % Generating a stochastic voltage signal with model ACF
    X = randn(nSamples,nPulseTrains);
    V = B*X;
    % Calculating the unbiased ACF of the stochastic voltage signal
    xCell = mat2cell(V,nSamples,ones(1,nPulseTrains));
    [RArray] = cellfun(@(x) xcorr(x,'unbiased'),xCell,'UniformOutput',false);
%     [RArray,lagNumber] = cellfun(@(x) xcorr(x,'unbiased'),xCell,'UniformOutput',false);
%     lagNumber = reshape(cell2mat(lagNumber),[nSamples*2-1,nPulseTrains]);
%     lagNumber(2:end+1,:)=lagNumber(1:end,:);
%     lagNumber(1,:)=lagNumber(2,:)-1;
    RArray = cell2mat(RArray);

    % Adding an ACF value at an additional to make the ACF length 2*nSamples
    % (instead of 2*nSamples-1)
    RArray(2:end+1,:) = RArray(1:end,:);
    RArray(1,:)=RArray(2,:);
    
    %Reorganizing ACF to match the input lag array
    RArray = ifftshift(RArray);
    RArray3D(:,:,thisExperiment)=RArray;
end

%% Step 3: Estimating the covariance and fitting the ACF
for thisPulseTrain = 1:1:nPulseTrains
    % Finding the covariance of the ACF by using samples from all
    % experiments
    covR = cov(squeeze(RArray3D(:,thisPulseTrain,:))');% Signal, PulseTrain, Experiment
    stdDeviation = (diag(covR)).^0.5;
    
%   data = mean(RArray,2);
    for thisExperiment = 1:1:nExperiments
        data=RArray3D(:,thisPulseTrain,thisExperiment);
        % Creating the function handle F(x): (model-data)./standard_deviation
        F = @(inputNLR) (autocorrelation_model...
            (inputNLR(1),inputNLR(2),inputNLR(4),inputNLR(3),lag)-data)./stdDeviation;
        %   Levenberg-Marquardt fitting
        outputTemp = lsqnonlin(F,inputNLR0,[],[],options);
        % Storing the fitted parameter values in the output variable
        outputNLR(:,thisPulseTrain,thisExperiment) = outputTemp; 

    end
end
%% Step:4 Statistical distribution of fitted parameters
% Generating the matrices that contain fitted parameters for each
% pulsetrain and experiment/ensemble member
    noiseDist = squeeze(outputNLR(1,:,:));
    signalDist = squeeze(outputNLR(2,:,:));
    dopplerShiftDist = squeeze(outputNLR(3,:,:));
    fWidthDist = squeeze(outputNLR(4,:,:));
     
end

% Defining the model ACF
function RmArray = autocorrelation_model(N,A,fWidth,omegaDoppler,lag)
       RmArray =N*double(lag==0) + A.*exp(-abs(lag).*fWidth).*exp(-1i*omegaDoppler.*lag);
   end

