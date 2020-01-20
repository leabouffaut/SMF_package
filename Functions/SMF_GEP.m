% h = SM_GEP(vecs, covs, covb, Qmax)
%
% Generation of the SMF filter bank based on the signal and noise 
% covariances matrices.
% The filters are the one maximizing the output SNR, wich is equivalent to
% solving the Generalized eigenvalue problem putting in relation the signal
% and noise covariance matrices
%
% INPUTS
%   - covs : Signal covariance matrix 
%   - vecs : Eigenvectors of the signal covariance matrix
%   - covb : Noise covariance matrix
%   - Qmax : Maximal number of filter (selected number of eigenvalues)
%
% OUTPUTS
%   - h : SNR-maximazing filter bank
% 
% Part of this code originally written by G. Julien, PhD (2010)


function h = SMF_GEP(vecs, covs, covb, Qmax)
% The signal and noise covariance matrices are projected onto a new
% subspace (the signal's subspace) and their order is reducted to Qmax (as
% allowed by the Karhunen Loeve theorem 
M = vecs(:,1:Qmax)'*covs*vecs(:,1:Qmax); % Projection of the signal covariance matrix
R = vecs(:,1:Qmax)'*covb*vecs(:,1:Qmax); % Projection of the noise covariance matrix

% Generalized eigenvalue problem
[vec,Lambda] = eig(M/R);
Lambda = abs(diag(Lambda));
[Lambda,pos] = sort(Lambda,'descend'); % Sort by descending lambda
vec = vec(:,pos); clear val

% Estimation of Phi and Psi
% The two matrices are orthogonal
% Phi is solution to the GEP
PHI = vecs(:,1:Qmax)*vec;
for n = 1:size(PHI,2),
    PHI(:,n) = PHI(:,n)/sqrt(PHI(:,n)'*covb*PHI(:,n)); % Normalisation
end
PSI = covb*PHI; 

% Filter bank
h = zeros(size(PHI));
Qmax = length(Lambda);
[N,~] = size(PHI);
for m = 2:Qmax+1
    h(:,m) = h(:,m-1) + PSI(:,m-1).*PHI(fix((N+1)/2),m-1);
end
h=h(:,2:m);
