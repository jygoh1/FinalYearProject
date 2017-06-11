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

fileNames = datasample(fileNames,floor(length(fileNames)/10),'Replace',false);

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

snrMDKFmask = zeros(numSNR, numFiles);
snrMaskLPC = zeros(numSNR, numFiles);
snrNoiseMask = zeros(numSNR, numFiles);
snrMDKF = zeros(numSNR, numFiles);
snrMMSE = zeros(numSNR, numFiles);
snrNoisy = zeros(numSNR, numFiles);

pesqMDKFmask = zeros(numSNR, numFiles);
pesqLPCmask = zeros(numSNR, numFiles);
pesqNoiseMask = zeros(numSNR, numFiles);
pesqMDKF = zeros(numSNR, numFiles);
pesqMMSE = zeros(numSNR, numFiles);
pesqNoisy = zeros(numSNR, numFiles);

stoiMDKFmask = zeros(numSNR, numFiles);
stoiLPCmask = zeros(numSNR, numFiles);
stoiNoiseMask = zeros(numSNR, numFiles);
stoiMDKF = zeros(numSNR, numFiles);
stoiMMSE = zeros(numSNR, numFiles);
stoiNoisy = zeros(numSNR, numFiles);

%%
for k = 1:numFiles
    % read in a speech file
    [clean,fs] = readsph(fileNames{k},'wt');

    clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

    ns = length(clean);       % number of speech samples
    v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

    for i = 1:length(targetSNR)
        noisy = v_addnoise(clean, fs, targetSNR(i), 'nzZ', v);  % add noise at chosen level keeping speech at 0 dB
        [snrMDKFmask(i,k), snrMaskLPC(i,k), snrNoiseMask(i,k), snrMDKF(i,k), snrMMSE(i,k), snrNoisy(i,k), ...
            pesqMDKFmask(i,k), pesqLPCmask(i,k), pesqNoiseMask(i,k), pesqMDKF(i,k), pesqMMSE(i,k), pesqNoisy(i,k), ...
            stoiMDKFmask(i,k), stoiLPCmask(i,k), stoiNoiseMask(i,k), stoiMDKF(i,k), stoiMMSE(i,k), stoiNoisy(i,k)] = runAll(clean, noisy, fs);
    end
    k
    
    save('snrStats_mmseLPCs','snrMDKFmask','snrMaskLPC','snrNoiseMask','snrMDKF','snrMMSE','snrNoisy');
    save('pesqStats_mmseLPCs','pesqMDKFmask','pesqLPCmask','pesqNoiseMask','pesqMDKF','pesqMMSE','pesqNoisy');
    save('stoiStats_mmseLPCs','stoiMDKFmask','stoiLPCmask','stoiNoiseMask','stoiMDKF','stoiMMSE','stoiNoisy');
end

%%
snrMDKFmask = snrMDKFmask(:,1:k);
snrMaskLPC = snrMaskLPC(:,1:k);
snrNoiseMask = snrNoiseMask(:,1:k);
snrMDKF = snrMDKF(:,1:k);
snrMMSE = snrMMSE(:,1:k);
snrNoisy = snrNoisy(:,1:k);

pesqMDKFmask = pesqMDKFmask(:,1:k);
pesqLPCmask = pesqLPCmask(:,1:k);
pesqNoiseMask = pesqNoiseMask(:,1:k);
pesqMDKF = pesqMDKF(:,1:k);
pesqMMSE = pesqMMSE(:,1:k);
pesqNoisy = pesqNoisy(:,1:k);

stoiMDKFmask = stoiMDKFmask(:,1:k);
stoiLPCmask = stoiLPCmask(:,1:k);
stoiNoiseMask = stoiNoiseMask(:,1:k);
stoiMDKF = stoiMDKF(:,1:k);
stoiMMSE = stoiMMSE(:,1:k);
stoiNoisy = stoiNoisy(:,1:k);

save('snrStats_mmseLPCs','snrMDKFmask','snrMaskLPC','snrNoiseMask','snrMDKF','snrMMSE','snrNoisy');
save('pesqStats_mmseLPCs','pesqMDKFmask','pesqLPCmask','pesqNoiseMask','pesqMDKF','pesqMMSE','pesqNoisy');
save('stoiStats_mmseLPCs','stoiMDKFmask','stoiLPCmask','stoiNoiseMask','stoiMDKF','stoiMMSE','stoiNoisy');
