function [S,f] = get_power_spectra_1(R,sampleTime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fs = 1/sampleTime;
L = length(R);
NFFT = 2^(nextpow2(L)+1);
% Power normalized by dividing by the frequency range
S = fftshift((2*real(fft(R,NFFT))-R(1,1)))*1/Fs;
f = Fs/2*(linspace(-1,1-2/NFFT,NFFT));

end

