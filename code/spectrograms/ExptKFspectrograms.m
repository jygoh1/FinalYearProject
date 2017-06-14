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

y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB

ns = length(y_clean);       % number of speech samples

noises = {'white'};
targetSNR = -20;      % if noiselevel = 5, target SNR is 5dB

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,targetSNR,'nzZ',v); % add noise at chosen level keeping speech at 0 dB

%%
load 'maskstats';

Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
LC = 0;             % LC for IBM (in dB)
p_MDKF = 2;
Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;
maskth = 0.15;
maskth_noise = 1.0;


mode='pJimc';
% mode='pJic';


%% spectrogram BMMDKF
m=2;
n=2;


figure;

subplot(m,n,1);

spgrambw(y_clean, fs, mode);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');
hold on
plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
hold off


subplot(m,n,2);

spgrambw(y_babble, fs, mode);
get(gca,'XTickLabel');
title(['Speech corrupted with white noise (' num2str(targetSNR) ' dB SNR)']);


subplot(m,n,3);

y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);
spgrambw(y_MDKF, fs, mode);
title(['MDKF with p=' num2str(p_MDKF)]);
hold on
plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
hold off

subplot(m,n,4);

y_BMMDKF = uncorrelatedMDKF_IBM_all(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, u_present, var_present, u_absent, var_absent);
spgrambw(y_BMMDKF, fs, mode);
title(['BMMDKF with p=' num2str(p_MDKF)]);
hold on
plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
hold off


% %% spectrogram LMDKF
% m=2;
% n=2;
% 
% figure;
% 
% subplot(m,n,1);
% 
% spgrambw(y_clean, fs, mode);
% get(gca,'XTickLabel');
% title('Spectrogram of clean speech');
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off
% 
% 
% subplot(m,n,2);
% 
% spgrambw(y_babble, fs, mode);
% get(gca,'XTickLabel');
% title(['Speech corrupted with white noise (' num2str(targetSNR) ' dB SNR)']);
% 
% 
% subplot(m,n,3);
% 
% y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);
% spgrambw(y_MDKF, fs, mode);
% title(['MDKF with p=' num2str(p_MDKF)]);
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off
% 
% 
% subplot(m,n,4);
% 
% y_LMDKF = MDKFmaskLPC(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, maskth);
% spgrambw(y_LMDKF, fs, mode);
% title(['LMDKF with p=' num2str(p_MDKF)]);
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off


% %% spectrogram NMDKF
% m=2;
% n=2;
% 
% figure;
% 
% subplot(m,n,1);
% 
% spgrambw(y_clean, fs, mode);
% get(gca,'XTickLabel');
% title('Spectrogram of clean speech');
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off
% 
% 
% subplot(m,n,2);
% 
% spgrambw(y_babble, fs, mode);
% get(gca,'XTickLabel');
% title(['Speech corrupted with white noise (' num2str(targetSNR) ' dB SNR)']);
% 
% 
% subplot(m,n,3);
% 
% y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);
% spgrambw(y_MDKF, fs, mode);
% title(['MDKF with p=' num2str(p_MDKF)]);
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off
% 
% 
% subplot(m,n,4);
% 
% y_NMDKF = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, maskth_noise);
% spgrambw(y_NMDKF, fs, mode);
% title(['NMDKF with p=' num2str(p_MDKF)]);
% hold on
% plot(get(gca,'xlim'), [1.5e3 1.5e3], 'b-.'); % Adapts to x limits of current axes
% hold off
