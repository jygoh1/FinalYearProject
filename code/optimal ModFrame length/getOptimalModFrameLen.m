close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));

%% generate mask statistics on set of training files

% acoustic frame
Tw = 16e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)

% MDKF order
p_MDKF = 2;

% modulation frame
Tw_slow = [12 16 20 24 32 48 64]*1e-3;
Ts_slow = 4e-3;
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
fileNames = datasample(fileNames,floor(length(fileNames)/8),'Replace',false);

%%
% assuming all fs same (which they are)
% then applying same windowing function gives same STFT matrix dimensions
% i.e. same number of freq bins

numLens = length(Tw_slow);      % number of modulation frame lengths to try
numFiles = length(fileNames);   % number of files to average over

[clean,fs] = readsph(fileNames{1},'wt');

excitationVarAllFiles = zeros(numLens,numFiles);

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

for i = 1:numLens
    
    [W_slow, ~, Ns_slow] = makeHammingWindow(fs_slow, Tw_slow(i), Ts_slow);
    
    for k = 1:numFiles
        % read in a speech file
        [clean,fs] = readsph(fileNames{k},'wt');
        clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

        excitationVarAllFiles(i,k) = getExcitationVar(clean, W, Nw, Ns, p_MDKF, W_slow, Ns_slow);
    end
end

excitationVar = nanmean(excitationVarAllFiles,2);

%%
len = [];
for i = 1:length(Tw_slow)
    len{i} = num2str(Tw_slow(i)*1e3);
end

bar(excitationVar)

set(gca, 'XTick', 1:length(Tw_slow), 'XTickLabel', len)
aaaa = get(gca,'XTickLabel');
title('\fontsize{23}Plot of excitation error vs. modulation frame lengths');
xlabel('\fontsize{22}Modulation Frame Length (ms)');
ylabel({'\fontsize{22}Average excitation error'});
set(gca,'XTickLabel',aaaa,'fontsize',20)