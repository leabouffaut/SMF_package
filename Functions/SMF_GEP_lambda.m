% Lambda = SMF_GEP_lambda(vecs, covs, covb, Qmax)
%
% Estimation of the eigenvalues of the Generalized eigenvalue problem putting 
% in relation the signal and online noise covariance matrices
%
% INPUTS
%   - covs : Signal covariance matrix 
%   - vecs : Eigenvectors of the signal covariance matrix
%   - covb : Noise covariance matrix
%   - Qmax : Maximal number of filter (selected number of eigenvalues)
%
% OUTPUTS
%   - Lambda : eigenvalues
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


function Lambda = SMF_GEP_lambda(vecs, covs, covb, Qmax)
% The signal and noise covariance matrices are projected onto a new
% subspace (the signal's subspace) and their order is reducted to Qmax (as
% allowed by the Karhunen Loeve theorem 
M = vecs(:,1:Qmax)'*covs*vecs(:,1:Qmax);
R = vecs(:,1:Qmax)'*covb*vecs(:,1:Qmax);

% Generalized eigenvalue problem
[~,Lambda] = eig(M/R);
Lambda = abs(diag(Lambda));
[Lambda,~] = sort(Lambda,'descend');