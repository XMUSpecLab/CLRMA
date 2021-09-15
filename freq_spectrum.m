function [f,p] = freq_spectrum(X,Fs)
% function of compute the frequency spectrum
% L表示信号长度
L = length(X);
% Fs表示采样频率固定设为2000,并非实际意义上的采样频率
if nargin == 1
    Fs = 2;
end
% 则信号周期为T
T = 1/Fs;
% 频率域的定义f
f = Fs*(0:floor(L/2))/L;
% 对信号进行FFT
Y = fft(X);
% 计算双边频谱
P = abs(Y/L);
% 计算单边频谱
p = P(1:floor(L/2)+1);
p(2:end-1) = 2*p(2:end-1);
end
