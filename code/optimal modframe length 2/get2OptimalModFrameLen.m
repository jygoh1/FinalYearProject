close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

warning('off','all')
warning

%% generate mask statistics on set of training files

% acoustic frame
Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)

% MDKF order
p_MDKF = 2;

% modulation frame
Tw_slow = [16 20 24 32 64 128 256]*1e-3;
Ts_slow = Tw_slow/4;
fs_slow = 1/Ts;

%% get excitation variance over range of modulation frame lengths
% averaged over multiple speech samples

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

%%
fileNames = datasample(fileNames,floor(length(fileNames)/12),'Replace',false);

%%
% assuming all fs same (which they are)
% then applying same windowing function gives same STFT matrix dimensions
% i.e. same number of freq bins

numLens = length(Tw_slow);      % number of modulation frame lengths to try
numFiles = length(fileNames)    % number of files to average over

[clean,fs] = readsph(fileNames{1},'wt');

excitationVar_other_AllFiles = zeros(numLens,numFiles);

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

pesqMDKF = zeros(numLens,numFiles);
stoiMDKF = zeros(numLens,numFiles);

noises = {'white'};
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);

for i = 1:numLens
    for k = 1:numFiles
        % read in a speech file
        [clean,fs] = readsph(fileNames{k},'wt');
        clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

        ns = length(clean);       % number of speech samples
        v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

        noisy = v_addnoise(clean, fs, 0, 'nzZ', v);     % white noise at 0 dB SNR for testing
        y_MDKF = idealMDKF_linear(noisy, clean, fs, Tw, Ts, p_MDKF, Tw_slow(i), Ts_slow(i), fs_slow);
        
        audiowrite('FYP\testfiles\clean.wav',clean,fs);
        audiowrite('FYP\testfiles\y_MDKF.wav',y_MDKF,fs);
        
        pesqMDKF(i,k) = pesqITU(fs,'FYP\testfiles\clean.wav','FYP\testfiles\y_MDKF.wav');

        cutoff = min([length(clean), length(y_MDKF)]);
        clean = clean(1:cutoff);
        y_MDKF = y_MDKF(1:cutoff);
        stoiMDKF(i,k) = stoi(clean,y_MDKF,fs);
    end
    i
    save('pesqMDKF','pesqMDKF');
    save('stoiMDKF','stoiMDKF');
end

pesqMDKF_avg = nanmean(pesqMDKF,2);
stoiMDKF_avg = nanmean(stoiMDKF,2);

%%
save('pesqMDKF_avg','pesqMDKF_avg');
save('stoiMDKF_avg','stoiMDKF_avg');

%%
load('pesqMDKF');
load('stoiMDKF');
load('pesqMDKF_avg');
load('stoiMDKF_avg');

len = [];
for i = 1:length(Tw_slow)
    len{i} = num2str(Tw_slow(i)*1e3);
end


figure;
plot(pesqMDKF_avg,'linewidth',1.5);

set(gca, 'XTick', 1:length(Tw_slow), 'XTickLabel', len)
title('\fontsize{19}Average PESQ values');
xlabel('\fontsize{14}Modulation Frame Length (ms)');
ylabel({'\fontsize{14}PESQ'});
aaaa = get(gca,'XTickLabel');
set(gca,'XTickLabel',aaaa,'fontsize',12)


figure;
plot(stoiMDKF_avg,'linewidth',1.5);

set(gca, 'XTick', 1:length(Tw_slow), 'XTickLabel', len)
title('\fontsize{19}Average STOI values');
xlabel('\fontsize{14}Modulation Frame Length (ms)');
ylabel({'\fontsize{14}STOI'});
aaaa = get(gca,'XTickLabel');
set(gca,'XTickLabel',aaaa,'fontsize',12)