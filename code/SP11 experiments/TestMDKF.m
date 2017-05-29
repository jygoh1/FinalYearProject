close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));

%% generate clean and noise signals

databases = '\\sapfs.ee.ic.ac.uk\Databases\';
timit = [databases 'Speech\TIMIT\TIMIT\TRAIN\'];
nato = [databases 'Noises\NatoNoise0\'];

% read in a speech file
[y_clean,fs,words,phonemes] = readsph([timit 'DR1\MDPK0\SX333.wav'],'wt');
% [y_clean,fs] = audioread('sp15.wav');

% downsample from 16 -> 8kHz
fs = fs/2;
y_clean = downsample(y_clean,2);

y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB

ns = length(y_clean);       % number of speech samples

noises = {'white'};
noiselevel = -5;      % if noiselevel = 5, target SNR is -5dB

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,-noiselevel,'nzZ',v); % add noise at chosen level keeping speech at 0 dB

%% spectrogram
m=2;
n=2;

mode='pJwiat';   % see spgrambw for list of modes
bw = [];
fmax = [];
dbrange = [-12 20];     % power in dB to plot on spectrogram
tinc = [];

figure;

subplot(m,n,1);
spgrambw(y_clean, fs, mode);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');


subplot(m,n,2);
spgrambw(y_babble, fs, mode);
get(gca,'XTickLabel');
title(['Speech corrupted with white Gaussian noise (' num2str(noiselevel) 'dB SNR)']);


subplot(m,n,4);
p_MDKF = 2;
Tw = 32e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
y_MDKF = idealMDKF(y_babble, y_clean, fs, Tw, Ts, p_MDKF);       % reduce noise but has distortion
spgrambw(y_MDKF, fs, mode);
title(['MDKF with p=' num2str(p_MDKF)]);