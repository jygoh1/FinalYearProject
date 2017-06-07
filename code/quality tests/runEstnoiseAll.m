function [pesqRatio, stoiRatio, segSNRratio] = runEstnoiseAll(y_clean, y_babble, fs, estnoisemask_th)

Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
LC = 0;             % LC for IBM (in dB)
p_MDKF = 2;
Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;

%% clean LPCs
% linear MDKF
y_MDKF = idealMDKF_linear(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);

% estnoise modified by IBM
y_MDKF_estnoisemask = idealMDKF_noiseIBM(y_babble, y_clean, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, estnoisemask_th);

%% PESQ
audiowrite('FYP\testfiles\y_clean.wav',y_clean,fs);
audiowrite('FYP\testfiles\y_MDKF.wav',y_MDKF,fs);
audiowrite('FYP\testfiles\y_MDKF_estnoisemask.wav',y_MDKF_estnoisemask,fs);

pesqMDKF = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF.wav');
pesqNoisemask = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_estnoisemask.wav');

pesqRatio = pesqNoisemask/pesqMDKF;

%% segSNR
cutoff = min([length(y_clean), length(y_MDKF), length(y_MDKF_estnoisemask)]);

y_clean = y_clean(1:cutoff);
y_MDKF = y_MDKF(1:cutoff);
y_MDKF_estnoisemask = y_MDKF_estnoisemask(1:cutoff);

segSNRMDKF = snrseg(y_MDKF, y_clean, fs);
segSNRMDKFmask = snrseg(y_MDKF_estnoisemask, y_clean, fs);

segSNRratio = segSNRMDKFmask/segSNRMDKF;

%% STOI
stoiMDKF = stoi(y_clean,y_MDKF,fs);
stoiNoisemask = stoi(y_clean,y_MDKF_estnoisemask,fs);

stoiRatio = stoiNoisemask/stoiMDKF;

end