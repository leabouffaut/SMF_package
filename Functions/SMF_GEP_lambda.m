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