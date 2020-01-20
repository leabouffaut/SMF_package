% [xcut, fs, N_begin_time, N_end_time] = cutfile_generalized(name, begin_time, duration)
%
% Load only the portion of file indicated by name that starts at begin_time 
% and lasts duration.
%
% INPUTS
%   - name : name and directory to the file 
%   - begin_time :  time of the begining of the selection (heure.min)
%          eg. 00.34  21.00  10.10   08.55
%   - duration : duration of the selection in MINUTES
%
% OUTPUTS
%   - xcut : selected signal vector
%   - fs : Sampling frequency (Hz)
%   - N_begin_time : number of samples between the first sample of the
%   audiofile and the begining of the selection
%   - N_end_time : number of samples between the first sample of the
%   audiofile and the end of the selection
% 
% WARNING: this works with Matlab 2015b or later (for earlier versions 
% audioread replace by wavread)

function [xcut, fs, N_begin_time, N_end_time] = ...
                        cutfile_generalized(name, begin_time, duration)

[x,fs] = audioread(name);

x = x';

% Begin time in samples
begin_time_sec = floor(begin_time)*3600 + ...
    (begin_time-floor(begin_time))*60*100; % Begin time in seconds
N_begin_time = round(begin_time_sec * fs +1 ); % Begin time in samples

% Duration in samples
duration_sec = duration*60; % Selection duration in seconds 
N_duration = duration_sec *fs +1; % Selection duration in samples 

% Last sample of the selection
N_end_time = N_begin_time + N_duration; 

% Cut
xcut = x(N_begin_time :N_end_time);
