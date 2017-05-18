% evaluate estnoiseg
close all;
clear all;

databases='\\sapfs.ee.ic.ac.uk\Databases\';
timit=[databases 'Speech\TIMIT\TIMIT\TRAIN\'];
nato=[databases 'Noises\NatoNoise0\'];
tinc=16e-3;         % frame hop = 16 ms
ovf=2;              % overlap factor = 2
% noises={'white', 'babble', 'factory1'};
noises={'white'};
nnoise=length(noises);  % number of noise signals
noiselev=-30:2:2;      % range of noise levels
nnlev=length(noiselev); % number of noise levels

% read in a speech file

[s,fs]=readsph([timit 'DR1\FCJF0\SA1.wav']);
% soundsc(s,fs) % optionally play the clean speech
s=activlev(s,fs,'n');   % normalize active level to 0 dB
ns=length(s);  % number of speech samples
ninc=round(0.016*fs);   % frame increment [fs=sample frequency]
win=sqrt(hanning(ovf*ninc,'periodic')); % define the window

% now read in the noises

v=zeros(ns,nnoise); % space for the noises
for j=1:nnoise
    [vj,fsj]=readwav([nato noises{j}]);
    vjr = resample(vj,fs,fsj);
    % soundsc(vir(1:3*fs),fs) % optionally play 3s of the resampled noise
    v(:,j)=vjr(1:ns)/std(vjr(1:ns)); % extract the initial chunck of noise and set to 0 dB
end

pow=zeros(nnlev,nnoise);     % space for the noises level estimates
for j=1:nnoise
for i=1:nnlev
    y=v_addnoise(s,fs,-noiselev(i),'nzZ',v(:,j)); % add noise at chosen level keeping speech at 0 dB
     % soundsc(y,fs) % optionally play 3s of the noisy speech
    f=enframe(y,win,ninc,'sp');
    x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
    pow(i,j)=10*log10(sum(mean(x,1),2));
end
end
figure;
plot(noiselev,pow,noiselev,noiselev,':k');
xlabel('True noise power (dB) [speech = 0 dB]');
ylabel('Estimated Noise Power (dB)');
legend(noises,'location','northwest');
figure;
plot(noiselev,pow-repmat(noiselev',1,nnoise),noiselev([1 end]),[0 0],':k');
xlabel('True noise power (dB) [speech = 0 dB]');
ylabel('Estimation error (dB)');
legend(noises,'location','northeast');


MODE='pJ';      % append 'c' to include colour bar
m=2;
n=2;

figure;

subplot(m,n,1);
spgrambw(s, fs, MODE);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y, fs, MODE);
get(gca,'XTickLabel');
title(['Spectrogram of speech embedded in multitalker babble at -' num2str(noiselev(i)) 'dB SNR']);
