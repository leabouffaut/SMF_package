% [TF,nu] = leafft(X,fs);
% Retourne la valeur absolue de la fft
% - X signal
% - fs Frequence ?chantillonage
% - TF r?sultat de la fft
% - nu fr?quences


function [TF,nu] = leafft(X,fs,n)
N = length(X);
if nargin == 2
    n = N;
else
    N = n;
end
TF = abs(fftshift(fft(X,n)))/N;
%if mod(N,2)==0
    dnu = fs/(N-1); 
%else
%dnu=fs/(N-1);
%end
nu = -fs/2:dnu:fs/2 ;





