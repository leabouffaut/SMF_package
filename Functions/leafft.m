% [TF,nu] = leafft(X,fs,n)
%
% Two sided abs(FFT) with frequency axis
% 
% INPUTS :
%   -  x: signal to analyze
%   - fs: sampling frequency (Hz)
%   -  n: (optional) FFT size 
% OUTPUTS :
%   - FTX: Two sided abs(FFT)
%   - nu: frequency axis


function [FTX,nu] = leafft(x,fs,n)
N = length(x);
if nargin == 2
    n = N;
else
    N = n;
end
FTX = abs(fftshift(fft(x,n)))/N;
dnu = fs/(N-1); 
nu = -fs/2:dnu:fs/2 ;





