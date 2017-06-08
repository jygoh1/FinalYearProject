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

fileNames = datasample(fileNames,floor(length(fileNames)/12),'Replace',false);

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

noisepow_orig = zeros(numSNR, numFiles);
noisepow_IBM = zeros(numSNR, numFiles);

Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
LC = 0;             % LC for IBM (in dB)
p_MDKF = 2;
Tw_slow = 24e-3;     % window and shift for each KF (in seconds)
Ts_slow = 4e-3;
fs_slow = 1/Ts;

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

%%
for k = 1:numFiles
    % read in a speech file
    [clean,fs] = readsph(fileNames{k},'wt');

    clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

    ns = length(clean);       % number of speech samples
    v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

    clean_fft = rfft(enframe(clean,W,Ns),Nw,2);
    
    for i = 1:length(targetSNR)
        noisy = v_addnoise(clean, fs, targetSNR(i), 'nzZ', v);  % add noise at chosen level keeping speech at 0 dB

        noisy_fft = rfft(enframe(noisy,W,Ns),Nw,2);     % STFT of noisy signal
        
        f_noisy = enframe(noisy,W,Ns,'sp');
        x_orig = estnoiseg(f_noisy,Ns/fs);   % estimate the noise power spectrum; each row is 1 time frame
        noisepow_orig(i,k) = 10*log10(sum(mean(x_orig,1),2));
        
        [~, mask] = IdBM(noisy, clean, fs, Tw, Ts, LC);
        x_IBM = estnoiseIBM(f_noisy,Ns/fs,mask',1.0);       % optimal mask threshold of 1.0
        noisepow_IBM(i,k) = 10*log10(sum(mean(x_IBM,1),2));
    end
    k
end

%%

save('estnoise_SNR','noisepow_orig');
save('estnoiseIBM_SNR','noisepow_IBM');
