% function [covs,vecs] = SMF_sig_preprocess(p)
%
% Estimation of the signal covariance matrix and eigenvectors
%
% INPUTS : 
%   - p : temporal reference signal (wavform)
%
% OUTPUTS :
%   - covs : signal covariance matrix
%   - vecs : signal sorted eigenvectors
%
% Piece of code originally written by G. Julien, PhD (2010)

function [covs,vecs] = SMF_sig_preprocess(p)

% Estimation of the signal covariance matrix covs
N = length(p);
autos = xcorr(p-mean(p)); % centered
autos = autos/max(autos); % between 1 and 0
covs = toeplitz(autos(N:end)); % To mirror diagonal

% Signal covariance matrix eigenvectors
[vecs,vp] = eig(covs); 
vp = diag(vp);
[~,pos] = sort(vp,'descend');
vecs = vecs(:,pos); % Sorted eigenvectors
