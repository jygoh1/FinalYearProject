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

y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);       % reduce noise but has distortion

y_NMDKF_1 = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th, 1.5);       % reduce noise but has distortion
y_NMDKF_2 = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th, 15);       % reduce noise but has distortion
y_NMDKF_3 = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th, 35);       % reduce noise but has distortion
y_NMDKF_4 = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th, 60);       % reduce noise but has distortion
y_NMDKF_5 = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th, 100);       % reduce noise but has distortion


%% PESQ
audiowrite('FYP\testfiles\y_clean.wav',y_clean,fs);
audiowrite('FYP\testfiles\y_MDKF.wav',y_MDKF,fs);
audiowrite('FYP\testfiles\y_NMDKF_1.wav',y_NMDKF_1,fs);
audiowrite('FYP\testfiles\y_NMDKF_2.wav',y_NMDKF_2,fs);
audiowrite('FYP\testfiles\y_NMDKF_3.wav',y_NMDKF_3,fs);
audiowrite('FYP\testfiles\y_NMDKF_4.wav',y_NMDKF_4,fs);
audiowrite('FYP\testfiles\y_NMDKF_5.wav',y_NMDKF_5,fs);

pesqMDKFnoiseIBM_1 = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_NMDKF_1.wav');
pesqMDKFnoiseIBM_2 = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_NMDKF_2.wav');
pesqMDKFnoiseIBM_3 = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_NMDKF_3.wav');
pesqMDKFnoiseIBM_4 = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_NMDKF_4.wav');
pesqMDKFnoiseIBM_5 = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_NMDKF_5.wav');
pesqMDKF = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF.wav');

fprintf(['PESQ ratio 1: ' num2str(pesqMDKFnoiseIBM_1/pesqMDKF) '\n']);
fprintf(['PESQ ratio 2: ' num2str(pesqMDKFnoiseIBM_2/pesqMDKF) '\n']);
fprintf(['PESQ ratio 3: ' num2str(pesqMDKFnoiseIBM_3/pesqMDKF) '\n']);
fprintf(['PESQ ratio 4: ' num2str(pesqMDKFnoiseIBM_4/pesqMDKF) '\n']);
fprintf(['PESQ ratio 5: ' num2str(pesqMDKFnoiseIBM_5/pesqMDKF) '\n']);

%% segSNR
cutoff = min([length(y_clean), length(y_MDKF), length(y_NMDKF_1)]);

y_clean = y_clean(1:cutoff);
y_MDKF = y_MDKF(1:cutoff);
y_NMDKF_1 = y_NMDKF_1(1:cutoff);
y_NMDKF_2 = y_NMDKF_2(1:cutoff);
y_NMDKF_3 = y_NMDKF_3(1:cutoff);
y_NMDKF_4 = y_NMDKF_4(1:cutoff);
y_NMDKF_5 = y_NMDKF_5(1:cutoff);

segSNRMDKFnoiseIBM_1 = snrseg(y_NMDKF_1, y_clean, fs);
segSNRMDKFnoiseIBM_2 = snrseg(y_NMDKF_2, y_clean, fs);
segSNRMDKFnoiseIBM_3 = snrseg(y_NMDKF_3, y_clean, fs);
segSNRMDKFnoiseIBM_4 = snrseg(y_NMDKF_4, y_clean, fs);
segSNRMDKFnoiseIBM_5 = snrseg(y_NMDKF_5, y_clean, fs);
segSNRMDKF = snrseg(y_MDKF, y_clean, fs);

fprintf(['segSNR ratio 1: ' num2str(segSNRMDKFnoiseIBM_1/segSNRMDKF) '\n']);
fprintf(['segSNR ratio 2: ' num2str(segSNRMDKFnoiseIBM_2/segSNRMDKF) '\n']);
fprintf(['segSNR ratio 3: ' num2str(segSNRMDKFnoiseIBM_3/segSNRMDKF) '\n']);
fprintf(['segSNR ratio 4: ' num2str(segSNRMDKFnoiseIBM_4/segSNRMDKF) '\n']);
fprintf(['segSNR ratio 5: ' num2str(segSNRMDKFnoiseIBM_5/segSNRMDKF) '\n']);

%% STOI
stoiMDKF = stoi(y_clean,y_MDKF,fs);
stoiMDKF_IBMnoise_1 = stoi(y_clean,y_NMDKF_1,fs);
stoiMDKF_IBMnoise_2 = stoi(y_clean,y_NMDKF_2,fs);
stoiMDKF_IBMnoise_3 = stoi(y_clean,y_NMDKF_3,fs);
stoiMDKF_IBMnoise_4 = stoi(y_clean,y_NMDKF_4,fs);
stoiMDKF_IBMnoise_5 = stoi(y_clean,y_NMDKF_5,fs);

fprintf(['STOI ratio 1: ' num2str(stoiMDKF_IBMnoise_1/stoiMDKF) '\n']);
fprintf(['STOI ratio 2: ' num2str(stoiMDKF_IBMnoise_2/stoiMDKF) '\n']);
fprintf(['STOI ratio 3: ' num2str(stoiMDKF_IBMnoise_3/stoiMDKF) '\n']);
fprintf(['STOI ratio 4: ' num2str(stoiMDKF_IBMnoise_4/stoiMDKF) '\n']);
fprintf(['STOI ratio 5: ' num2str(stoiMDKF_IBMnoise_5/stoiMDKF) '\n']);
