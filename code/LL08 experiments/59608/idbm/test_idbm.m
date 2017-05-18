% Ideal binary mask test framework by Kamil Wojcicki, 2011 (test_idbm.m)
clear all; close all; clc; randn('seed',0); rand('seed',0); fprintf('.\n');


    % specify clean speech input wav file
    file.clean = 'sp10.wav';

    % read clean speech samples from wav file 
    [ speech.clean, fs, nbits ] = wavread(file.clean);

    % create time vector
    time = [ 0:length(speech.clean)-1 ]/fs; 

    % read noise samples from wav file 
    % [ noise, fs, nbits ] = wavread( 'ssn.wav' );
    [ noise, fs, nbits ] = wavread( 'babble.wav' );

    SNR = 0; % noisy speech signal-to-noise ratio (dB)

    % mix clean speech and noise at a prescribed SNR
    speech.noisy = addnoise( speech.clean, noise, SNR ); 

    Tw = 32;    % analysis frame duration (ms) 
    Ts = Tw/8;  % analysis frame shift (ms)
    LC = -5;    % local SNR criterion (LC) 

    % apply ideal binary mask to the noisy speech signal 
    speech.processed = idbm( speech.noisy, speech.clean, fs, Tw, Ts, LC ); 

    % get stimuli labels and their count 
    methods = fieldnames(speech);
    M = length(methods);

    % loop through stimuli and compute SegSNR scores
    system(sprintf('rm -f ./%s.txt', mfilename));
    diary(sprintf('%s.txt', mfilename)); diary on;
    for m = 1:M 
        method = methods{m};
        mos.(method) = segsnr( speech.clean, speech.(method), fs );
        fprintf('SegSNR %12s : %5.2f\n', method, mos.(method));
    end
    diary off;

    % plot stimuli
    figure('Position', [20 20 1000 200*M], 'PaperPositionMode', 'auto', 'Visible', 'on');
    for m = 1:M 
        method = methods{m};

        % time domain plots
        subplot(M,2,2*m-1); 
        plot(time,speech.(method),'k-'); 
        xlim([min(time) max(time)]);
        title(sprintf('Waveform: %s  SegSNR=%0.2f', method, mos.(method)), 'interpreter', 'none');
        xlabel('Time (s)');
        ylabel('Amplitude');

        % spectrogram plots
        subplot(M,2,2*m); 
        myspectrogram(speech.(method), fs);
        set(gca,'ytick',[0:1000:16000],'yticklabel',[0:16]);
        title(sprintf('Spectrogram: %s  SegSNR=%0.2f', method, mos.(method)), 'interpreter', 'none');
        xlabel('Time (s)');
        ylabel('Frequency (kHz)');
    end
    print('-dpng', sprintf('%s.png', mfilename));

    % write stimuli samples to wav files
    for m = 1:M 
        method = methods{m};
        audio.(method) = 0.999*speech.(method)./max(abs(speech.(method)));
        wavwrite(audio.(method), fs, nbits, sprintf('%s.wav',method));
    end


% EOF
