close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

warning('off','all')
warning

%%
databases = '\\sapfs.ee.ic.ac.uk\Databases\';
timit = [databases 'Speech\TIMIT\TIMIT\TRAIN\'];
db = [timit 'DR1\'];
nato = [databases 'Noises\NatoNoise0\'];

folders = dir(db);

for k = length(folders):-1:1
    % remove non-folders
    if ~folders(k).isdir
        folders(k) = [ ];
        continue
    end

    % remove folders starting with .
    fname = folders(k).name;
    if fname(1) == '.'
        folders(k) = [ ];
    end
end

fileNames = {};
for i=1:length(folders)
    folder = folders(i).name;
    directory = [db folder];
    files = dir(fullfile(directory,'*.wav'));
    files = {files.name}';                      %'# file names

    data = cell(numel(files),1);                %# store file contents
    for j=1:numel(files)
        fname = fullfile(directory,files{j});     %# full path to file
        fileNames = [fileNames, fname];
    end
end

fileNames = datasample(fileNames,floor(length(fileNames)/15),'Replace',false);

%%
[clean,fs] = readsph(fileNames{1},'wt');

numFiles = length(fileNames)
noises = {'white'};

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);

targetSNR = [-20 -15 -10 -5 0 5 10 15 20];
targetSNR = targetSNR';

numSNR = length(targetSNR);

pesqMaskLPC_minuspt20 = zeros(numSNR, numFiles);
pesqMaskLPC_minuspt15 = zeros(numSNR, numFiles);
pesqMaskLPC_minuspt10 = zeros(numSNR, numFiles);
pesqMaskLPC_minuspt05 = zeros(numSNR, numFiles);
pesqMaskLPC_pt05 = zeros(numSNR, numFiles);
pesqMaskLPC_pt10 = zeros(numSNR, numFiles);
pesqMaskLPC_pt15 = zeros(numSNR, numFiles);
pesqMaskLPC_pt20 = zeros(numSNR, numFiles);
pesqMaskLPC_pt25 = zeros(numSNR, numFiles);

stoiMaskLPC_minuspt20 = zeros(numSNR, numFiles);
stoiMaskLPC_minuspt15 = zeros(numSNR, numFiles);
stoiMaskLPC_minuspt10 = zeros(numSNR, numFiles);
stoiMaskLPC_minuspt05 = zeros(numSNR, numFiles);
stoiMaskLPC_pt05 = zeros(numSNR, numFiles);
stoiMaskLPC_pt10 = zeros(numSNR, numFiles);
stoiMaskLPC_pt15 = zeros(numSNR, numFiles);
stoiMaskLPC_pt20 = zeros(numSNR, numFiles);
stoiMaskLPC_pt25 = zeros(numSNR, numFiles);

Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
LC = 0;             % LC for IBM (in dB)
p_MDKF = 2;
Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;

%%
for k = 1:numFiles
    [clean,fs] = readsph(fileNames{k},'wt');

    clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

    ns = length(clean);       % number of speech samples
    v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

    for i = 1:length(targetSNR)
        noisy = v_addnoise(clean, fs, targetSNR(i), 'nzZ', v);  % add noise at chosen level keeping speech at 0 dB
        
        y_mmse = ssubmmse(noisy, fs);

        % linear MDKF
        y_MDKF = idealMDKF_linear(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow);

        % mask-enhanced LPCs
        y_MDKF_lpcMask_minuspt20 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, -0.2);
        y_MDKF_lpcMask_minuspt15 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, -0.15);
        y_MDKF_lpcMask_minuspt10 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, -0.1);
        y_MDKF_lpcMask_minuspt05 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, -0.05);
        y_MDKF_lpcMask_pt05 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, 0.05);
        y_MDKF_lpcMask_pt10 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, 0.1);
        y_MDKF_lpcMask_pt15 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, 0.15);
        y_MDKF_lpcMask_pt20 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, 0.2);
        y_MDKF_lpcMask_pt25 = MDKFmaskLPC(noisy, y_mmse, fs, Tw, Ts, p_MDKF, Tw_slow, Ts_slow, fs_slow, LC, 0.25);
        
        audiowrite('FYP\testfiles\y_clean.wav',y_clean,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_minuspt20.wav',y_MDKF_lpcMask_minuspt20,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_minuspt15.wav',y_MDKF_lpcMask_minuspt15,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_minuspt10.wav',y_MDKF_lpcMask_minuspt10,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_minuspt05.wav',y_MDKF_lpcMask_minuspt05,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_pt05.wav',y_MDKF_lpcMask_pt05,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_pt10.wav',y_MDKF_lpcMask_pt10,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_pt15.wav',y_MDKF_lpcMask_pt15,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_pt20.wav',y_MDKF_lpcMask_pt20,fs);
        audiowrite('FYP\testfiles\y_MDKF_lpcMask_pt25.wav',y_MDKF_lpcMask_pt25,fs);

        pesqMaskLPC_minuspt20(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_minuspt20.wav');
        pesqMaskLPC_minuspt15(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_minuspt15.wav');
        pesqMaskLPC_minuspt10(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_minuspt10.wav');
        pesqMaskLPC_minuspt05(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_minuspt05.wav');
        pesqMaskLPC_pt05(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_pt05.wav');
        pesqMaskLPC_pt10(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_pt10.wav');
        pesqMaskLPC_pt15(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_pt15.wav');
        pesqMaskLPC_pt20(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_pt20.wav');
        pesqMaskLPC_pt25(i,k) = pesqITU(fs,'FYP\testfiles\y_clean.wav','FYP\testfiles\y_MDKF_lpcMask_pt25.wav');

        pesqMaskLPC_minuspt20(i,k) = stoi(y_clean,y_MDKF_lpcMask_minuspt20,fs);
        pesqMaskLPC_minuspt15(i,k) = stoi(y_clean,y_MDKF_lpcMask_minuspt15,fs);
        pesqMaskLPC_minuspt10(i,k) = stoi(y_clean,y_MDKF_lpcMask_minuspt10,fs);
        pesqMaskLPC_minuspt05(i,k) = stoi(y_clean,y_MDKF_lpcMask_minuspt05,fs);
        pesqMaskLPC_pt05(i,k) = stoi(y_clean,y_MDKF_lpcMask_pt05,fs);
        pesqMaskLPC_pt10(i,k) = stoi(y_clean,y_MDKF_lpcMask_pt10,fs);
        pesqMaskLPC_pt15(i,k) = stoi(y_clean,y_MDKF_lpcMask_pt15,fs);
        pesqMaskLPC_pt20(i,k) = stoi(y_clean,y_MDKF_lpcMask_pt20,fs);
        pesqMaskLPC_pt25(i,k) = stoi(y_clean,y_MDKF_lpcMask_pt25,fs);
        
    end
    k
end
