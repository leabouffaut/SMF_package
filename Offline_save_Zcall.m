%% Simulate Z-call sigmoid and estimate signal's covariance matrix
%
% In this program:
% (1) First, the Z-call sigmoid is simulated based on FX. Socheleau 2015 
% JASA paper.
% (2) then, the covariance matrix of the signal is evaluated as a first 
% part of the SMF application
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


clearvars
close all
clc

addpath Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Save = 'On';
disp(['Save is: ' Save])
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sigmoid parameters
fs = 100 ; % (Hz) sampling frequency

Tz = 20 ; % (s) duree d'un zcall
tz = 0:1/fs:Tz; % time axis
N = length(tz); % number of samples

fc = 22.6 ;  % (Hz) Z-call center frequency
alpha = 1.8; % slope

% Following parameters, 3 values are considered to account for Z-call
% frequency variations
L = [-4.5 -4 -3.5]; % Lower assymptote
U = [3.2 3.6 4];    % Upper asymptote 
M = [Tz/2 (Tz+0.5)/2 (Tz+1)/2];  % (s) time at half of the slope


% %% Instantaneous frequency
% This section is not mandatory, it can just be used to visualize the
% result simulated call
%
% f_whale = zeros(length(L),N);
% nb_samples_tshift = floor(tshift*fs);
% for i = 1:length(L)
%     f_whale(i,:) = fc + L(i) + ((U(i)-L(i))./(1+exp(alpha*(tz-M(i)))));
% end

% Estimation of the time-varying phase 
adj = fc - 8.5; % Compared to the Socheleau paper, the sigmoid has to be 
                % shifted in frequency to match the Z-call signal
L = L-adj;
U = U-adj;

phase_whale = zeros(1,N);
s_whale = zeros(1,N);
i = 1;
for n = 0:N-1
    phase_whale(n+1) = 2*pi*(L(i)*n/fs + ((U(i)-L(i))/alpha)*log((1+exp(-alpha*M(i)))/(1+exp(alpha*(n/fs-M(i)))))); % Expression de laphase
end
phase_whale = phase_whale(end:-1:1);
s_whale = s_whale + exp(1j*phase_whale); % temporal signal creation


% Amplitude
s_whale = s_whale/max(abs(s_whale)); 
amplitude = [ones(1,round(N/2)-1) ones(1,round(N/2))*0.95];
p0 = 1e-6;
s_whale = amplitude.*s_whale*p0;
s_whale = real(s_whale/max(abs(real(s_whale))));


%% Calcul de la martice de covariance du signal
[covs,vecs] = SMF_sig_preprocess(real(s_whale)) ;

if strcmp(Save,'On') == 1
save('Offline_saved/s_whale.mat', 's_whale','covs','vecs')
disp('Matrices saved in: Offline_saved/s_whale.mat')
else disp('Matrices not saved')
end


