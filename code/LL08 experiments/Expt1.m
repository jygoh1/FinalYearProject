close all;
clearvars;
clc

addpath(genpath('FYP'));
addpath(genpath('voicebox'));

%% generate clean and noise signals

numFiles=1;

y_clean=[];
for i=1:numFiles
    clean = sprintf('sp0%d.wav',i);
    [y_next,fs] = audioread(clean);
    y_clean = [y_clean;y_next];
end

noiseSNR=5;
y_babble=[];
for i=1:numFiles
    babble = sprintf('sp0%d_babble_sn%d.wav',i,noiseSNR);
    [y_next,fs] = audioread(babble);
    y_babble = [y_babble;y_next];
end

% y_babble=v_addnoise(y_clean,fs,noiseSNR);      % this seems to add much more noise than necessary, results are MUCH worse than in paper

% sound(y_clean,fs);
% sound(y_babble,fs);

%% spectrogram

MODE='pJ';      % append 'c' to include colour bar
m=4;
n=2;

figure;
subplot(m,n,1);
spgrambw(y_clean, fs, MODE);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y_clean, fs, MODE);
get(gca,'XTickLabel');
title('Spectrogram of clean speech');

subplot(m,n,3);
spgrambw(y_babble, fs, MODE);
get(gca,'XTickLabel');
title(['Spectrogram of speech embedded in multitalker babble at -' num2str(noiseSNR) 'dB SNR']);

subplot(m,n,4);
spgrambw(y_babble, fs, MODE);
get(gca,'XTickLabel');
title(['Spectrogram of speech embedded in multitalker babble at -' num2str(noiseSNR) 'dB SNR']);

subplot(m,n,5);
Tw=20;          % frame duration in ms
overlap=50/100; % 50% overlap
Ts=Tw*overlap;          % frame shift in ms (overlap)
LC=-10;         % local SNR criterion in dB
[y_ibmOutput, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC );       % reduce noise but has distortion
[T,F,B]=spgrambw(y_ibmOutput', fs, MODE);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title(['IBM with SNR threshold of ' num2str(LC) 'dB']);

subplot(m,n,7);
spgrambw(y_ibmOutput, fs, MODE);
get(gca,'XTickLabel');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) 'dB']);

subplot(m,n,6);
LC=0;         % local SNR criterion in dB
[y_ibmOutput, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC );       % reduce noise but has distortion
[T,F,B]=spgrambw(y_ibmOutput', fs, MODE);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title(['IBM with SNR threshold of ' num2str(LC) 'dB']);

subplot(m,n,8);
spgrambw(y_ibmOutput, fs, MODE);
get(gca,'XTickLabel');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) 'dB']);

% sound(y_ibmOutput,fs);

%% spectrogram

MODE='pJ';      % append 'c' to include colour bar
m=4;
n=1;

figure;
subplot(m,n,1);
spgrambw(y_clean, fs, MODE);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title('Spectrogram of clean speech');

subplot(m,n,2);
spgrambw(y_babble, fs, MODE);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of speech embedded in multitalker babble at -' num2str(noiseSNR) 'dB SNR']);

subplot(m,n,3);
LC=0;         % local SNR criterion in dB
[y_ibmOutput, mask] = IdBM(y_babble, y_clean, fs, Tw, Ts, LC );       % reduce noise but has distortion
[T,F,B]=spgrambw(y_ibmOutput', fs, MODE);
imagesc(T,F/1e3,mask);         % plot it manually
axis('xy');
get(gca,'XTickLabel');
xlabel('Time (s)');
ylabel('Freq. (kHz)');
title(['IBM with SNR threshold of ' num2str(LC) 'dB']);

subplot(m,n,4);
spgrambw(y_ibmOutput, fs, MODE);
get(gca,'XTickLabel');
ylabel('Freq. (kHz)');
title(['Spectrogram of segregated mixture using IBM SNR threshold of ' num2str(LC) 'dB']);

%% intelligibility test
x_axis = [-50,-40,-30,-20,-10,0,5,10,15];
fiveDbResults = [25,40,80,92,95,95,91,77,35];
tenDbResults = [10,15,70,90,94,95,80,40,12];

figure;
plot(x_axis,fiveDbResults,'o-')
hold on;
plot(x_axis,tenDbResults,'*-')
hold on;
plot([0 0],[1,100],'k--')
% yL = get(gca,'YLim');
% line([0 0],yL,'Color','k-');

labels = {'UN','-40','-30','-20','-10','0','10','20'};
set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
title('Speech Intelligibility Test of IBM');
legend({'-5dB SNR','-10dB SNR'});
xlabel('LC Threshold (dB)');
ylabel({'Percentage of words identified correctly'});