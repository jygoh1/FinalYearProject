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

y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB

% zeropad = fs*0.5;
% y_clean = padarray(y_clean,[zeropad 0]);      % zero-extend by approx 500ms for reliable noise estimation

ns = length(y_clean);       % number of speech samples

noises = {'white'};
targetSNR = 5;      % if noiselevel = 5, target SNR is -5dB

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,targetSNR,'nzZ',v); % add noise at chosen level keeping speech at 0 dB

%% spectrogram
m=5;
n=2;

mode='pJcwiatl';      % append 'c' to include colour bar
bw = [];
fmax = [];
dbrange = [-12 20];     % power in dB to plot on spectrogram
tinc = [];

figure;
subplot(m,n,1);
spgrambw(y_clean, fs, mode, bw, fmax, dbrange, tinc, words);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y_clean, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title('Spectrogram of clean speech');

subplot(m,n,3);
spgrambw(y_babble, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of speech embedded in multitalker babble at ' num2str(targetSNR) ' dB SNR']);

subplot(m,n,4);
spgrambw(y_babble, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of speech embedded in multitalker babble at ' num2str(targetSNR) ' dB SNR']);

subplot(m,n,5);
Tw = 20e-3;                 % frame duration in s
overlap = 50/100;           % 50% overlap
Ts = Tw*overlap;            % frame shift in s (overlap)
LC = -10;                   % local SNR criterion in dB
[y_IBM_1, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC);       % reduce noise but has distortion
[T,F,B] = spgrambw(y_IBM_1', fs, mode, bw, fmax, dbrange, tinc, phonemes);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
xlabel('Time (s)');
ylabel('Freq. (kHz)');
title(['IBM with SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,7);
spgrambw(y_IBM_1, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,9);
[y_TBM_1, TBMmask] = TBM(y_babble, y_clean, fs, Tw, Ts, LC ); 
spgrambw(y_TBM_1, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using TBM SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,6);
LC = 0;         % local SNR criterion in dB
[y_IBM_2, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC );       % reduce noise but has distortion
[T,F,B] = spgrambw(y_IBM_2', fs, mode, bw, fmax, dbrange, tinc, phonemes);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
xlabel('Time (s)');
ylabel('Freq. (kHz)');
title(['IBM with SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,8);
spgrambw(y_IBM_2, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,10);
[y_TBM_2, TBMmask] = TBM(y_babble, y_clean, fs, Tw, Ts, LC ); 
spgrambw(y_TBM_2, fs, mode, bw, fmax, dbrange, tinc, phonemes);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using TBM SNR threshold of ' num2str(LC) ' dB']);


%% spectrogram
m=4;
n=1;

mode='pJim';      % append 'c' to include colour bar

figure;
subplot(m,n,1);
spgrambw(y_clean, fs, mode);
get(gca,'XTickLabel');
xlabel('');
% ylabel('Freq. (kHz)');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y_babble, fs, mode);
get(gca,'XTickLabel');
xlabel('');
% ylabel('Freq. (kHz)');
title(['Speech corrupted with white Gaussian noise at ' num2str(targetSNR) ' dB SNR']);

subplot(m,n,3);
LC = 0;         % local SNR criterion in dB
[y_IBM, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC );       % reduce noise but has distortion
[T,F,B] = spgrambw(y_IBM', fs, mode);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
ylabel('Frequency (kMel)');
title(['IBM with SNR threshold of ' num2str(LC) ' dB']);

subplot(m,n,4);
spgrambw(y_IBM, fs, mode);
get(gca,'XTickLabel');
% ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) ' dB']);