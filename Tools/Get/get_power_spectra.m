function [S,f,moments] = get_power_spectra(x,t,mode)
%UNTITLED2 Summary of this function goes here
%   x - 'signal':signal or 'acf':ACF 
%   t - 'signal':time   or 'acf':time 
%   note t>0
%   mode - 'signal'(default) Inputs are signal, time
%        - 'acf' Inputs are positive-lag ACF, lag, and dt
%   dt - sample time
if nargin < 3
    mode =0;
end
    
if strcmp(mode,'signal')
    dt = t(2)-t(1); % Sample Time
    nSamples = length(t);
%     nEnsembles=size(x,2);
%     xCell = mat2cell(x,nSamples,ones(1,nEnsembles));
%     [RArray,lagNumber] = cellfun(@xcorr,xCell,'UniformOutput',false);
%     RArray = cell2mat(RArray);
%     lagNumber = reshape(cell2mat(lagNumber),[nSamples*2-1,nEnsembles]);
%     RUnbiasedCoeff = xcorr(ones(nSamples,1));
%     RUnbiasedCoeff = repmat(RUnbiasedCoeff,1,nEnsembles);
%     RArray = RArray./RUnbiasedCoeff;
%     lagNumber(end+1,:)=lagNumber(end,:)+1;
%     RArray(end+1,:)=RArray(end);
%     
%     
    
    [RArray,lagNumber] = xcorr(x,'unbiased');
    lagNumber(end+1)=lagNumber(end)+1; %To bring Rx to power of 2 samples
    RArray(end+1)=RArray(end);
        
    begLagIndx = find(lagNumber==0);
    [tempVal,endLagIndx] = max(lagNumber);
    endLagIndx=endLagIndx-1;
    lagIndx = begLagIndx:endLagIndx;

    RArrayDisplay = RArray(lagIndx);
    lagNumberDisplay = lagNumber(lagIndx);
    lagDisplay = lagNumberDisplay.*dt;

    RArrayModified = ifftshift(RArray);
    lagNumberShifted = ifftshift(lagNumber); %X

elseif strcmp(mode,'acf')
    dt = t(2)-t(1); % Sample Time
    nSamples = length(t);
    RArrayDisplay = x;
    lagDisplay = t;
    
    endIndx = nSamples;
    RArrayModified = zeros(1,2*endIndx);
    RArrayModified(1,1:endIndx) = x(1:end);
    RArrayModified(1,endIndx+1) = RArrayModified(endIndx);
    RArrayModified(1,endIndx+2:end) = flipud(conj(x(2:end))); 
    RArrayModified = RArrayModified';
end    
    
    
R = toeplitz(conj(RArrayDisplay));  %Traditional Autocorrelation matrix
R1 = toeplitz(conj(RArrayModified)); %Contains reflection of the RArray

% Calculate D, covR, covS.
D = dftmtx(length(RArrayModified)); %Question? Should the input be modified periodic ACF or just one-half

%PSD using the DFT matrix
S1 = fftshift(abs(D*RArrayModified)).*dt; %X

%PSD using the fft function
[S,f] = get_fft(RArrayDisplay,dt);

%Sample plot - X
% figure;
% plot(f/1000,S);
% hold on;
% plot(f/1000,S1,'.');

% Storing data
moments.ACF_Array=RArrayDisplay;
moments.ACF_Array_modified = RArrayModified;
moments.lag=lagDisplay;
moments.ACF=R;
moments.S_DFT = S1; %PSD from DFT matrix multiplication method
moments.D=D;

end

