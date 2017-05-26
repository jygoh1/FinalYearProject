% close all;
% clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));


%% generate clean and noise signals

databases = '\\sapfs.ee.ic.ac.uk\Databases\';
timit = [databases 'Speech\TIMIT\TIMIT\TRAIN\'];
nato = [databases 'Noises\NatoNoise0\'];

% read in a speech file
% [y_clean,fs] = readsph([timit 'DR1\FCJF0\SA1.wav']);
[y_clean,fs] = audioread('sp10.wav');
y_clean = activlev(y_clean,fs,'n');     % normalise active level to 0 dB


% % test with stationary signal (sine wave)
% t = [0:1/30000:0.1]';
% A = 1;
% fs = 6e3;
% y_clean = A*sin(2*pi*fs*t);


zeropad = fs*0.4;
y_clean = padarray(y_clean,[zeropad 0]);      % zero-extend by approx 500ms for reliable noise estimation

ns = length(y_clean);       % number of speech samples

noises = {'white'};
% noises = {'f16'};
noiselevel = -5;      % if noiselevel = 5, target SNR is -5dB

% read in the noises
[vj,fsj] = readwav([nato noises{1}]);
vjr = resample(vj,fs,fsj);
v = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB; v is noise

y_babble = v_addnoise(y_clean,fs,-noiselevel,'nzZ',v); % add noise at chosen level keeping speech at 0 dB


% y_babble = y_clean + randn(length(y_clean),1)*10^(noiselevel/10);


% for j=1:nnoise
%     [vj,fsj] = readwav([nato noises{j}]);
%     vjr = resample(vj,fs,fsj);
% %     soundsc(vjr(1:3*fs),fs) % optionally play 3s of the resampled noise
%     v(:,j) = vjr(1:ns)/std(vjr(1:ns));  % extract the initial chunck of noise and set to 0 dB
% end
% 
% pow=zeros(nnlev,nnoise);     % space for the noises level estimates
% for j=1:nnoise
%     for i=1:nnlev
%         y_babble = v_addnoise(y_clean,fs,-noiselevel(i),'nzZ',v(:,j)); % add noise at chosen level keeping speech at 0 dB
% %         soundsc(y_babble,fs) % optionally play 3s of the noisy speech
% %         f=enframe(y,win,ninc,'sp');
% %         x=estnoiseg(f,ninc/fs); % estimate the noise power spectrum
% %         pow(i,j)=10*log10(sum(mean(x,1),2));
%     end
% end


% figure;
% plot(noiselevel,pow,noiselevel,noiselevel,':k');
% xlabel('True noise power (dB) [speech = 0 dB]');
% ylabel('Estimated Noise Power (dB)');
% legend(noises,'location','northwest');
% figure;
% plot(noiselevel,pow-repmat(noiselevel',1,nnoise),noiselevel([1 end]),[0 0],':k');
% xlabel('True noise power (dB) [speech = 0 dB]');
% ylabel('Estimation error (dB)');
% legend(noises,'location','northeast');


%% generate clean and noise signals

% [y_clean,fs] = audioread('sp15.wav');
% 
% zeropad=fs*0.5;
% 
% noiselevel=-5;
% % y_clean=padarray(y_clean,[zeropad 0]);      % zero-extend by approx 500ms for reliable noise estimation
% % y_babble=v_addnoise(y_clean,fs,noiselevel);   % add white noise
% y_babble=y_clean+randn(length(y_clean),1)*10^(noiselevel/10);
% 
% 
% % soundsc(y_clean,fs);    % scale audio data to play as loudly as possible without clipping
% % sound(y_babble,fs);


%% spectrogram

mode='pJcwi';   % see spgrambw for list of modes
m=2;
n=2;

dbrange = [-10 20];     % power in dB to plot on spectrogram
bw = 200;   % default

figure;
subplot(m,n,1);
spgrambw(y_clean, fs, mode, bw, fs/2, dbrange);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y_babble, fs, mode, bw, fs/2, dbrange);
get(gca,'XTickLabel');
title(['Spectrogram of speech embedded in multitalker babble at -' num2str(noiselevel) 'dB SNR']);


p_TDKF=10;
% q_TDKF=4;

subplot(m,n,3);
Tw = 32e-3;         % frame duration in s
Ts = 4e-3;          % frame shift in s (overlap)
y_TDKF = TDKF(y_babble, y_clean, fs, Tw, Ts, p_TDKF);
% [T,F,B]=spgrambw(y_TDKF', fs, MODE);
% imagesc(T,F/1e3,B);         % plot it manually
% axis('xy');
% get(gca,'XTickLabel');
% xlabel('Time (s)');
% ylabel('Frequency (kHz)');
spgrambw(y_TDKF, fs, mode, bw, fs/2, dbrange);
title(['TDKF with p=' num2str(p_TDKF)]);


p_MDKF=2;

subplot(m,n,4);
Tw = 32e-3;         % frame duration in s  
Ts = 4e-3;          % frame shift in s (overlap)
y_MDKF = MDKF(y_babble, y_clean, fs, Tw, Ts, p_MDKF);       % reduce noise but has distortion
spgrambw(y_MDKF, fs, mode, bw, fs/2, dbrange);
title(['MDKF with p=' num2str(p_MDKF)]);


%%

% figure;
% subplot(m,n,1)
% plot(y_babble,'g')
% hold on;
% plot(y_TDKF,'r');
% hold on;
% plot(y_clean,'k')
% legend({'noisy signal','TDKF signal','Clean signal'});
% ylim([-8 8])
% 
% subplot(m,n,2)
% plot(y_TDKF,'k');
% title('TDKF-processed time-domain signal');
% ylim([-8 8])
% 
% subplot(m,n,3)
% plot(y_clean,'b')
% title('Clean time-domain signal');
% ylim([-8 8])
% 
% subplot(m,n,4)
% plot(y_babble,'b')
% title('noisy time-domain signal');
% ylim([-8 8])