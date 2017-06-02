close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

%%
%GETMASKSTATS Summary of this function goes here
%   Detailed explanation goes here

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

[clean,fs] = readsph(fileNames{1},'wt');

numFiles = length(fileNames)
noises = {'white'};

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);

targetSNR = [-20 -15 -10 -5 0 5 10 15 20];
targetSNR = targetSNR';

% numAlgos = 5;       % including noisy
numSNR = length(targetSNR);

snrMDKFmask = zeros(numSNR, numFiles);
snrMaskLPC = zeros(numSNR, numFiles);
snrMDKF = zeros(numSNR, numFiles);
snrMMSE = zeros(numSNR, numFiles);
snrNoisy = zeros(numSNR, numFiles);
pesqMDKF = zeros(numSNR, numFiles);
pesqMDKFmask = zeros(numSNR, numFiles);
pesqLPCmask = zeros(numSNR, numFiles);
pesqMMSE = zeros(numSNR, numFiles);
pesqNoisy = zeros(numSNR, numFiles);

for k = 1:numFiles
    % read in a speech file
    [clean,fs] = readsph(fileNames{k},'wt');

    clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

    ns = length(clean);       % number of speech samples
    v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

    for i = 1:length(targetSNR)
        noisy = v_addnoise(clean, fs, targetSNR(i), 'nzZ', v);  % add noise at chosen level keeping speech at 0 dB
        [snrMDKFmask(i,k), snrMaskLPC(i,k), snrMDKF(i,k), snrMMSE(i,k), snrNoisy(i,k), pesqMDKF(i,k), pesqMDKFmask(i,k), pesqLPCmask(i,k), pesqMMSE(i,k), pesqNoisy(i,k)] = runAll(clean, noisy, fs);
    end
    k
end

%%
snrMDKFmask = snrMDKFmask(:,1:k);
snrMaskLPC = snrMaskLPC(:,1:k);
snrMDKF = snrMDKF(:,1:k);
snrMMSE = snrMMSE(:,1:k);
snrNoisy = snrNoisy(:,1:k);
pesqMDKF = pesqMDKF(:,1:k);
pesqMDKFmask = pesqMDKFmask(:,1:k);
pesqLPCmask = pesqLPCmask(:,1:k);
pesqMMSE = pesqMMSE(:,1:k);
pesqNoisy = pesqNoisy(:,1:k);
save('snrStats_1','snrMDKFmask','snrMaskLPC','snrMDKF','snrMMSE','snrNoisy');
save('pesqStats_1','pesqMDKF','pesqMDKFmask','pesqLPCmask','pesqMMSE','pesqNoisy');

%%
load 'snrStats_1';
load 'pesqStats_1';


%%
% %% generate new test clean and noise signal
% load 'maskstats';
% 
% databases = '\\sapfs.ee.ic.ac.uk\Databases\';
% timit = [databases 'Speech\TIMIT\TIMIT\TEST\'];
% nato = [databases 'Noises\NatoNoise0\'];
% 
% % read in a speech file
% [y_clean,fs,words,phonemes] = readsph([timit 'DR1\MSJS1\SA1.wav'],'wt');
% 
% y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB
% 
% ns = length(y_clean);       % number of speech samples
% 
% noises = {'white'};
% 
% % read in the noises
% [vj,fsj] = readwav([nato noises{1}]);
% vjr = resample(vj,fs,fsj);
% v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise
% 
% %%
% targetSNR = [-5 0 5 10 20];
% targetSNR = targetSNR';
% 
% % numAlgos = 8;       % including noisy
% numAlgos = 5;       % including noisy
% segSNR = zeros(length(targetSNR), numAlgos);
% 
% for i = 1:length(targetSNR)
%     y_babble = v_addnoise(y_clean, fs, targetSNR(i), 'nzZ', v);  % add noise at chosen level keeping speech at 0 dB
%     segSNR(i,:) = runAll(y_clean, y_babble, fs);
% end
% 
% 
% %%
% close all;
% 
% figure;
% h = plot(targetSNR,segSNR(:,[1,3,4,5]),'linewidth',1.2);
% 
% title('\fontsize{23}Average segSNR values');
% xlabel('\fontsize{17}Global SNR of noisy speech (dB)');
% ylabel({'\fontsize{17}segSNR (dB)'});
% aaaa = get(gca,'XTickLabel');
% set(h,{'Marker'},{'x';'s';'o';'*'})
% set(gca,'XTickLabel',aaaa,'fontsize',15)
% legendCell = {'MDKF-IBM', 'MDKF', 'MMSE', 'Noisy'};
% legend(legendCell,'FontSize',12)
% 
% figure;
% h = plot(targetSNR,segSNR(:,[2,3,4,5]),'linewidth',1.2);
% 
% title('\fontsize{23}Average segSNR values');
% xlabel('\fontsize{17}Global SNR of noisy speech (dB)');
% ylabel({'\fontsize{17}segSNR (dB)'});
% aaaa = get(gca,'XTickLabel');
% set(h,{'Marker'},{'s';'o';'*';'^'})
% set(gca,'XTickLabel',aaaa,'fontsize',15)
% % legendCell = {'MDKF', 'MDKF_uncorr', 'MDKF_IBM', 'MDKF_uncorrIBM', 'LPC-enhanced', 'MMSE', 'Noisy'};
% legendCell = {'LPC-enhanced', 'MDKF', 'MMSE', 'Noisy'};
% legend(legendCell,'FontSize',12)
% 
% % %%
% % figure;
% % h = plot(targetSNR,segSNR(:,[6,8]),'linewidth',1.2);
% % 
% % title('\fontsize{23}Average segSNR values');
% % xlabel('\fontsize{17}Global SNR of noisy speech (dB)');
% % ylabel({'\fontsize{17}segSNR (dB)'});
% % aaaa = get(gca,'XTickLabel');
% % set(h,{'Marker'},{'s';'^'})
% % set(gca,'XTickLabel',aaaa,'fontsize',15)
% % legendCell = {'MDKF-IBM (clean)', 'MDKF (clean)'};
% % legend(legendCell,'FontSize',12)
% % 
% % figure;
% % h = plot(targetSNR,segSNR(:,[7,8]),'linewidth',1.2);
% % 
% % title('\fontsize{23}Average segSNR values');
% % xlabel('\fontsize{17}Global SNR of noisy speech (dB)');
% % ylabel({'\fontsize{17}segSNR (dB)'});
% % aaaa = get(gca,'XTickLabel');
% % set(h,{'Marker'},{'s';'^'})
% % set(gca,'XTickLabel',aaaa,'fontsize',15)
% % % legendCell = {'MDKF', 'MDKF_uncorr', 'MDKF_IBM', 'MDKF_uncorrIBM', 'LPC-enhanced', 'MMSE', 'Noisy'};
% % legendCell = {'LPC-enhanced (clean)', 'MDKF (clean)'};
% % legend(legendCell,'FontSize',12)