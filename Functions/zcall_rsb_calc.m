% snr = zcall_rsb_calc(stft, stft_med, f, t, tx, N)
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

function snr = zcall_rsb_calc(stft,stft_med,f,t,tx,N)
% Step (1) Estimation of the presency of a Z-call
    flow = 25.8 ;
    fhigh = 26.6 ;
    a = find(f >= flow,1);
    b = find(f >= fhigh,1);

    snr_zcall = abs(max(stft(a:b,:))./mean(stft_med(a:b,:),1));
    [~,index_detect] = max(stft(a:b,:));
    freq_high = f(a+index_detect); freq_high = interp1(t,freq_high,tx);
    clear a b index
    snr_zcall = 20*log10(interp1(t,snr_zcall,tx));
    snr_zcall = medfilt1(snr_zcall,100);

% Step (2) estimation of the presency of a larger banc transient signal
    snr_ztrans = abs(stft./stft_med);
    a = find(f >= 15,1);
    b = find(f >= 27,1);
    snr_ztrans(a:b,:) = zeros(size(snr_ztrans(a:b,:)));
    snr_ztrans = max(snr_ztrans,[],1); clear a b
    snr_ztrans = 20*log10(interp1(t,snr_ztrans,tx));
    snr_ztrans = medfilt1(snr_ztrans,200);
    
% Step (3) SNR estimation
    snr = 20*log10((10.^(snr_zcall/20)./10.^(snr_ztrans/20)));
    if mean(snr_zcall(N:length(snr_zcall)-N)) - mean(snr_ztrans(N:length(snr_zcall)-N)) < 0
        seuil = 0;
    else
        seuil = mean(snr_zcall(N:length(snr_zcall)-N)) - mean(snr_ztrans(N:length(snr_zcall)-N));
    end
    
    snr = snr-seuil ;