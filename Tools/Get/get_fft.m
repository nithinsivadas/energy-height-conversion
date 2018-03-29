function [S,f] = get_fft(x,sampleTime,mode)
%get_fft: FFT of x(t) with t>=0 to generate S(f)
%   t - is the signal with t>=0
%   sampleTime - is the sampling time of the signal/ACF [s]
%   modes : 0 - (Default) Uses the FFT trick to calculate the sepctra
%                         only with the real part of the signal
%         : 1 - Reflects the signal to form a periodic signal before ffting
%         : 2 - Assumes the signal is already periodic
if nargin<3
    mode = 0;
end
Fs = 1/sampleTime;
L = length(x);
    if mode==0
        NFFT = 2^(nextpow2(L)+1);
        % Power normalized by dividing by the frequency range
        S = fftshift((2*real(fft(x,NFFT))-x(1,1)))*1/Fs;
        f = Fs/2*(linspace(-1,1-2/NFFT,NFFT));
    elseif mode==1
        endIndx = L;
        RArray = zeros(1,2*endIndx);
        RArray(1,1:endIndx) = x(1:end);
        RArray(1,endIndx+1) = RArray(endIndx);
        RArray(1,endIndx+2:end) = flipud(conj(x(2:end))); 
        NFFT = 2^nextpow2(length(RArray));
        S = abs(fftshift(fft(RArray,NFFT)))*1/Fs;
        f = Fs/2*(linspace(-1,1-2/NFFT,NFFT));
    elseif mode==2
        NFFT = 2^(nextpow2(L));
        S = abs(fftshift(fft(x,NFFT)))*1/Fs;
        f = Fs/2*(linspace(-1,1-2/NFFT,NFFT));
    else
        error('Mode does not match with the available options.');
    end
end

