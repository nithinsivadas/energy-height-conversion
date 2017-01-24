Fs = 3;            % Sampling frequency
T = 1/Fs;             % Sampling period
t = (1206518400.0000000:T:1206536400.0000000);        % Time vector
L = length(t);             % Length of signal

S = 2*sin(2*pi*(t-t(1)));
X=S;

figure;
subplot(2,1,1)
plot(1000*(t(1:50)-t(1)),X(1:50))
title('Signal')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y = abs(fft(X));

P2 = (abs(Y/L));
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
%subplot(2,1,2)
%plot(f,Y(1:L/2+1))

subplot(2,1,2)
plot(f,P1.^2);
%title('Single-Sided Amplitude Spectrum of X(t)')
%xlabel('f (Hz)')
%ylabel('|P1(f)|')