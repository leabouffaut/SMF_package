% [covn, Qmax,snr] = SMF_noise_rsb_preprocess_median_phase(x, fs, fft_size,overlap, N,med_win_size)
%
% Estimation of the ONLINE noise covariance matrix, Qmax and time-varying
% input SNR based onan estimation of the background noise;
% This estimation is realized by using the outlier-smoothing property of 
% the median filter, applying it to each frequency canal of the absolute
% STFT of the observation. In the STFT representation, the signal can be
% considered transient (occupy only a limited number of frequency bins) and
% can therefore be "smoothen out" by an adequately calibered median filter.
%
% INPUTS :
%   - x: observation where the background noise is estimated
%   - fs: sampling frequency (Hz) 
%   - fft_size: fft size for the spectrogram
%   - overlap: Percentage of window overlap for the spectrogram (%)
%   - N: length of the reference signal in bins 
%
% SORTIES :
%   - covn : Online-estimated noise covariance matrix
%   - Qmax : NMax number of filter
%   - snr : time-varying snr (dB)



function [covn, Qmax,snr] = SMF_noise_rsb_preprocess_median_phase(x, fs, fft_size,overlap, N,med_win_size)
% Time parameters
Tx = (length(x)-1)/fs; %(s) duration of the observation
tx = 0:1/fs:Tx; % temporal axis

% STFT calclation
[stft,f,t,~] = spectrogram(x,hann(fft_size),round((overlap/100)*fft_size),fft_size,fs);
stft_phase = angle(stft); % Keep the phase
stft_med = abs(stft); % absolute STFT to apply the medain filter
[~,Nt] = size(stft);

% Median-filtering
for j = 1 :Nt - med_win_size
    stft_med (:, j+floor(med_win_size/2)) = median(stft_med(:,j:j+med_win_size),2);
end


% SNR estimation
% WARNING: In this function, the SNR is estimated for the Z-call frequency
% band, if applied to an other signal, change frequency boundary inside the
% function zcall_rsb_calc
snr = zcall_rsb_calc(stft,stft_med,f,t,tx,N);

% estimation of the noise covariance matrix
autonoise = real(ifft(stft_med.^2,[],2)); 
autonoise = sum(autonoise,2);
autophase = sum(stft_phase.^2,2) ;                     
autonoise = autonoise/max(autonoise);
fnew = linspace(0,50,N);
autonoise = interp1(f,autonoise + 1j*autophase,fnew);  
covn = abs(toeplitz(autonoise)); % Cov du bruit   
covn = covn-min(min(covn));             
covn = covn/max(max(covn));             

% Analysis of the noise covariance matrix
[~,vp] = eig(covn); % Vecteurs propres vec
                    % Valeurs  propres vp
vp = diag(vp);
vp = sort(vp);
Qmax = sum(vp > 1);
