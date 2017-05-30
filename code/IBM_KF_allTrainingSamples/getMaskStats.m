function [u_present, var_present, u_absent, var_absent] = getMaskStats(noiselevel, Tw, Ts, LC)
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

%%
% assuming all fs same (which they are)
% then applying same windowing function gives same STFT matrix dimensions
% i.e. same number of freq bins

[clean,fs] = readsph(fileNames{1},'wt');
[W, Nw, Ns] = makeHammingWindow(fs/2, Tw, Ts);
test = rfft(enframe(clean,W,Ns),Nw,2);
numBins = size(test, 2);

numFiles = length(fileNames);
u_present_all = zeros(numBins,numFiles);
var_present_all = zeros(numBins,numFiles);
u_absent_all = zeros(numBins,numFiles);
var_absent_all = zeros(numBins,numFiles);

noises = {'white'};
    
for k=1:numFiles
    % read in a speech file
    [clean,fs] = readsph(fileNames{k},'wt');

    % downsample from 16 -> 8kHz
    fs = fs/2;
    clean = downsample(clean,2);

    clean = activlev(clean,fs,'n');     % normalise active level to 0 dB

    ns = length(clean);       % number of speech samples

    % read in the noises
    [vj,fsj] = readwav([nato noises{1}]);
    vjr = resample(vj,fs,fsj);
    v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

    noisy = v_addnoise(clean,fs,-noiselevel,'nzZ',v); % add noise at chosen level keeping speech at 0 dB

    [W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

    noisy_fft = rfft(enframe(noisy,W,Ns),Nw,2);     % STFT of noisy signal
    noisyfft_mag = abs(noisy_fft);
    noisyfft_mag = noisyfft_mag';
    
    % mask used to scale mean and variance of observation
    [~, mask] = IdBM(noisy, clean, fs, Tw, Ts, LC);
    for i = 1:numBins
        present = noisyfft_mag(i, mask(i,:) > 0);
        absent = noisyfft_mag(i, mask(i,:) == 0);
        u_present_all(i,k) = mean(present);
        var_present_all(i,k) = var(present);
        u_absent_all(i,k) = mean(absent);
        var_absent_all(i,k) = var(absent);
    end
end

% take average over large set of training samples, excluding NaN values
u_present = nanmean(u_present_all,2);
var_present = nanmean(var_present_all,2);
u_absent = nanmean(u_absent_all,2);
var_absent = nanmean(u_present_all,2);

end
