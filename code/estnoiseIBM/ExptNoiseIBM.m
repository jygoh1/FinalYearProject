close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));
warning('off','all')

%% generate clean and noise signals

databases = '\\sapfs.ee.ic.ac.uk\Databases\';
timit = [databases 'Speech\TIMIT\TIMIT\TRAIN\'];
nato = [databases 'Noises\NatoNoise0\'];

% read in a speech file
[y_clean,fs,words,phonemes] = readsph([timit 'DR1\MDPK0\SX333.wav'],'wt');
y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB

ns = length(y_clean);       % number of speech samples

noises = {'white'};
targetSNR = 0;      % if noiselevel = 5, target SNR is 5dB

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,targetSNR,'nzZ',v); % add noise at chosen level keeping speech at 0 dB

%% spectrogram
m=2;
n=2;

p_MDKF = 2;
Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)

Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;
LC = 0;

mask_th = 1;
mode='pJim'; 

% figure;
% 
% subplot(m,n,1);
% 
% spgrambw(y_clean, fs, mode);
% get(gca,'XTickLabel');
% title('Spectrogram of clean speech');
% 
% 
% subplot(m,n,2);
% 
% spgrambw(y_babble, fs, mode);
% get(gca,'XTickLabel');
% title(['Speech corrupted with white Gaussian noise (' num2str(targetSNR) ' dB SNR)']);


y_mmse = ssubmmse(y_babble, fs);

% subplot(m,n,3);

y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);       % reduce noise but has distortion
% spgrambw(y_MDKF, fs, mode);
% title(['MDKF with p=' num2str(p_MDKF)]);


% subplot(m,n,4);

y_MDKF_IBMnoise = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th);       % reduce noise but has distortion
% spgrambw(y_MDKF_IBMnoise, fs, mode);
% title(['IBM-noise MDKF with p=' num2str(p_MDKF)]);


%% PESQ
audiowrite('FYP\testfiles\y_clean.wav',y_clean,fs);
audiowrite('FYP\testfiles\y_MDKF.wav',y_MDKF,fs);
audiowrite('FYP\testfiles\y_MDKF_IBMnoise.wav',y_MDKF_IBMnoise,fs);

pesqMDKFnoiseIBM = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_IBMnoise.wav');
pesqMDKF = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF.wav');

pesqMDKFnoiseIBM/pesqMDKF

%% segSNR
cutoff = min([length(y_clean), length(y_MDKF), length(y_MDKF_IBMnoise)]);

y_clean = y_clean(1:cutoff);
y_MDKF = y_MDKF(1:cutoff);
y_MDKF_IBMnoise = y_MDKF_IBMnoise(1:cutoff);

segSNRMDKFnoiseIBM = snrseg(y_MDKF_IBMnoise, y_clean, fs);
segSNRMDKF = snrseg(y_MDKF, y_clean, fs);

segSNRMDKFnoiseIBM/segSNRMDKF

%% STOI
stoiMDKF = stoi(y_clean,y_MDKF,fs);
stoiMDKF_IBMnoise = stoi(y_clean,y_MDKF_IBMnoise,fs);

stoiMDKF_IBMnoise/stoiMDKF