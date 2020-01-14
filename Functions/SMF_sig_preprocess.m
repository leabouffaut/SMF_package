% function [covs,vecs] = SMF_sig_preprocess(p)
%
% 
% Full description of the method's theory is described in
%
%       L. Bouffaut, R. Dreo, V. Labat, A. Boudraa and G. Barruol 
%       'Passive stochastic matched filter for antarctic blue whale call
%       detection,' in J. Acoust. Soc. Am, 144(2) (2018).
%
% and any use of this material should refer accordingly. Part of the code
% was written by Gregory Julien, PhD (2010).
%
% Estimation of the signal covariance matrix and eigenvectors
%
% INPUTS : 
%   - p : temporal reference signal
%
% OUTPUTS :
%   - covs : signal covariance matrix
%   - vecs : signal sorted eigenvectors
%
% ------------------------------------------------------------------------%
% Author: Lea Bouffaut, PhD
% Naval Academy Research Institute, Brest, France.
% Date: 10-30-2019
% 
% Full description of the method's theory is described in
%
%       L. Bouffaut, R. Dreo, V. Labat, A. Boudraa and G. Barruol 
%       'Passive stochastic matched filter for antarctic blue whale call
%       detection,' in J. Acoust. Soc. Am, 144(2) (2018).
%
% and any use of this material should refer accordingly. Part of the code
% was written by Gregory Julien, PhD (2010).
% ------------------------------------------------------------------------%


function [covs,vecs] = SMF_sig_preprocess(p)

% Estimation of the signal covariance matrix covs
N = length(p);
autos = xcorr(p-mean(p)); % centered
autos = autos/max(autos); % between 1 and 0
covs = toeplitz(autos(N:end)); % To mirror diagonal

% Signal covariance matrix eigenvectors
[vecs,vp] = eig(covs);
vp = diag(vp);
[val,pos] = sort(vp,'descend');clear val
vecs = vecs(:,pos); % Sorted eigenvectors
