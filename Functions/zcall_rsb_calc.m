% [rsb,seuil,freq_high] = zcall_rsb_calc(stft,stft_med)
%
% Online SNR estimation
%
% INPUTS :
%  - stft : observation STFT
%  - stft_med : median-filtered STFT
%  - f : STFT frequency vector
%  - t :  STFT time vector
%  - tx : observation temporal vector for interpolation
%  - N : length of the reference signal in number of bins 
%
% OUTPUTS :
%  - snr : time varying snr (dB)
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
% and any use of this material should refer accordingly.
% ------------------------------------------------------------------------%

function snr = zcall_rsb_calc(stft,stft_med,f,t,tx,N)
% Step (1) Estimation of the presency of a Z-call
    fbasse = 25.8 ;
    fhaute = 26.6 ;
    a = find(f >= fbasse,1); % 26
    b = find(f >= fhaute,1); % 28

    rsb_zcall = abs(max(stft(a:b,:))./mean(stft_med(a:b,:),1));
    [~,index_detect] = max(stft(a:b,:));
    freq_high = f(a+index_detect); freq_high = interp1(t,freq_high,tx);
    clear a b index
    rsb_zcall = 20*log10(interp1(t,rsb_zcall,tx));
    rsb_zcall = medfilt1(rsb_zcall,100);

% Step (2) estimation of the presency of a larger banc transient signal
    rsb_ztrans = abs(stft./stft_med);
    a = find(f >= 15,1); %15
    b = find(f >= 27,1); %30
    rsb_ztrans(a:b,:) = zeros(size(rsb_ztrans(a:b,:)));
    rsb_ztrans = max(rsb_ztrans,[],1); clear a b
    rsb_ztrans = 20*log10(interp1(t,rsb_ztrans,tx));
    rsb_ztrans = medfilt1(rsb_ztrans,200);
    
% Step (3) SNR estimation
    snr = 20*log10((10.^(rsb_zcall/20)./10.^(rsb_ztrans/20)));
    if mean(rsb_zcall(N:length(rsb_zcall)-N)) - mean(rsb_ztrans(N:length(rsb_zcall)-N)) < 0
        seuil = 0;
    else
        seuil = mean(rsb_zcall(N:length(rsb_zcall)-N)) - mean(rsb_ztrans(N:length(rsb_zcall)-N));
    end
    
    snr = snr-seuil ;