close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

%%
Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
LC = 0;             % LC for IBM (in dB)
targetSNR = 0;      % -noise level (in dB)

%% generate new test clean and noise signal
load 'maskstats';

databases = '\\sapfs.ee.ic.ac.uk\Databases\';
timit = [databases 'Speech\TIMIT\TIMIT\TEST\'];
nato = [databases 'Noises\NatoNoise0\'];

% read in a speech file
[y_clean,fs,words,phonemes] = readsph([timit 'DR1\MSJS1\SA1.wav'],'wt');

y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB

ns = length(y_clean);       % number of speech samples

noises = {'white'};

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,targetSNR,'nzZ',v); % add noise at chosen level keeping speech at 0 dB


%% results
p_MDKF = 2;
Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;

mask_th = 0.5;      % for estimating LPCs

y_mmse = ssubmmse(y_babble, fs);

y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);

y_MDKF_uncorrIBM = uncorrelatedMDKF_IBM_all(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, u_present, var_present, u_absent, var_absent);

y_maskLPC = MDKFmaskLPC(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, mask_th);

%%
cutoff = min([length(y_clean), length(y_maskLPC), length(y_MDKF_uncorrIBM), length(y_MDKF)]);

y_clean = y_clean(1:cutoff);
y_MDKF = y_MDKF(1:cutoff);
y_MDKF_uncorrIBM = y_MDKF_uncorrIBM(1:cutoff);
y_maskLPC = y_maskLPC(1:cutoff);

mode = [];
segsnr_MDKF = snrseg(y_MDKF,y_clean,fs,mode);
segsnr_MDKFmask = snrseg(y_MDKF_uncorrIBM,y_clean,fs,mode);
segsnr_LPCmask = snrseg(y_maskLPC,y_clean,fs,mode);

%%
audiowrite('FYP\testfiles\y_clean.wav',y_clean,fs);
audiowrite('FYP\testfiles\y_MDKF.wav',y_MDKF,fs);
audiowrite('FYP\testfiles\y_MDKF_uncorrIBM.wav',y_MDKF_uncorrIBM,fs);
audiowrite('FYP\testfiles\y_maskLPC.wav',y_maskLPC,fs);

pesq_MDKF = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF.wav');
pesq_MDKFmask = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_uncorrIBM.wav');
pesq_LPCmask = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_maskLPC.wav');

pesq_LPCmask/pesq_MDKF
pesq_MDKFmask/pesq_MDKF

%%
stoi_MDKF = stoi(y_clean,y_MDKF,fs);
stoi_LPCmask = stoi(y_clean,y_maskLPC,fs);
stoi_MDKFmask = stoi(y_clean,y_MDKF_uncorrIBM,fs);