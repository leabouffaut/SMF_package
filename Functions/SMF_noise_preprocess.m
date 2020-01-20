% [covn, Qmax] = SMF_noise_preprocess(noise, N)
%
% Estimation of the noise covariance matrix an dof maximum number of 
% filters Qmax that can be applied.
%
% INPUTS :
%   - noise : noise waveform
%   - N : size (length in bins) of the reference signal
%
% OUTPUTS :
%   - covn : noise covariance matrix
%   - Qmax : maximum number of filters (max number of selected eigenvalues)
%
% Piece of code originally written by G. Julien, PhD (2010)

function [covn, Qmax] = SMF_noise_preprocess(noise, N)

% Noise covariance matrix estimation
autob = zeros(1,2*N-1);
for n=1:length(noise)-N+1,
    mini = noise(n:n+N-1); 
    provi = xcorr(mini-mean(mini)); % centered
    autob = autob+provi/max(provi); % reduced
end
autob = autob/max(autob); % reduced
covn = toeplitz(autob(N:end)); % Noise covariance matrix

% Noise covariance matrix analysis to evaluate Qmax
% Based on the highest eigenvalues
[~,valn] = eig(covn); 
valn = diag(valn);
valn = sort(valn);
Qmax = sum(valn > valn(fix(0.9*length(valn))));
% Keep 10% of the highest noise eigenvalues
% Noise highest eignevalues are representative of most of the noise energy 
