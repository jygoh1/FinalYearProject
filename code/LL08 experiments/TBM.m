function [y, mask] = TBM(noisy, clean, fs, Tw, Ts, LC)

% This function modifies the spectrum of a noisy input speech signal by 
% applying a target binary mask to the signal, assuming additive noise 
% corruption. The ideal binary mask is computed from an oracle (true) 
% signal-to-noise ratio (SNR) by thresholding with local SNR criterion
% specified in LC. The synthesized enhanced speech is returned in y
% and the binary mask is returned in mask.
% 99% for ideal, 75% for non-ideal
% 
% Inputs
%     noisy:    noisy speech signal (vector)
%     clean:    clean speech signal (vector)
%     fs:       sampling frequency (Hz)
%     Tw:       frame duration (s)
%     Ts:       frame shift (s)
%     LC:       local SNR criterion (dB)
% 
% Outputs
%     y:        enhanced speech signal i.e., noisy signal with the IBM applied
%     mask:     IBM
% 
% Author: Jia Ying Goh, Imperial College London, January 2017
% Credits:
%     Kamil Wojcicki, UTD, October 2011
%     Mike Brookes, Imperial College London - VOICEBOX: A speech processing toolbox for MATLAB

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);
F_clean = rfft(enframe(clean,W,Ns),Nw,2);      % do STFT: one row per time frame, +ve frequencies only
F_noisy = rfft(enframe(noisy,W,Ns),Nw,2);

% compute the true STFT noise spectrum (assumes additive noise distortion)
noise = F_clean - F_noisy;

% compute true SNR and threshold it to produce the ideal binary mask
% TBM formed by comparing energy in each T-F unit to percentage of total energy in signal
targetE = abs(F_clean).^2;
totalE = sum(sum(targetE)) / numel(targetE);
SNR = targetE / totalE;
% SNR = abs(F_clean).^2 ./ abs(noise).^2;     % IBM(t,f) = 1 if Target(t,f) - Masker(t,f) > LC (dB)
mask = zeros( size(SNR) );
mask( SNR>0.5*10^(0.1*LC) ) = 1;        % set values to 1 if SNR of signal higher than threshold SNR

% apply the ideal binary mask and create modified complex spectrum
F_clean_TBM = abs(F_noisy) .* mask .* exp(1i*angle(F_noisy));

y = overlapadd(irfft(F_clean_TBM,Nw,2),W,Ns);  % reconstitute the time waveform
mask = mask';

% EOF
