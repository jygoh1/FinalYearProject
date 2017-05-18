function [y, MASK] = IdBM( x, s, fs, Tw, Ts, LC )
% IDBM Apply spectral ideal binary mask to noisy speech.
%
%   [Y,MASK,SNR]=IDBM(X,S,FS,TW,TS,LC) performs FFT based short-time
%   spectral analysis-modification-synthesis, with frame duration 
%   of TW ms and frame shift of TS ms. Spectral modification 
%   performed is the application of an ideal binary mask to noise
%   corrupted speech signal, assuming additive noise distortion.
%   The ideal binary mask is computed from an oracle (true) signal-to-noise ratio (SNR) by
%   thresholding with local SNR criterion specified in LC. 
%   The synthesized enhanced speech is returned in Y, the ideal 
%   binary mask is returned in MASK, while the true instantaneous
%   spectral signal-to-noise ratio is stored in SNR.
%   
%   Inputs
%           X is the noisy speech signal as vector.
%
%           S is the clean speech signal as vector.
%
%           FS is the sampling frequency in Hz.
%
%           TW is the frame duration in ms.
%
%           TS is the frame shift in ms.
%           
%           LC is the local SNR criterion (dB).
%
%   Outputs 
%           Y is the enhanced speech signal, 
%           i.e., noisy signal with the ideal binary mask applied.
%
%           MASK is the ideal binary mask.
%
%   Author: Jia Ying Goh, Imperial College London, January 2017
%   Credits:
%       Kamil Wojcicki, UTD, October 2011
%       Mike Brookes, Imperial College London - VOICEBOX: A speech processing toolbox for MATLAB

    
    % check for correct number of input arguments
    if( nargin~=6 ), error( sprintf('Not enough input arguments. Type "help %s" for usage help.', mfilename) ); end;
    
    Nw = round( fs*Tw*0.001 );              % frame duration (in samples)
    Ns = round( fs*Ts*0.001 );              % frame shift (in samples)
    if Nw/Ns==4
        W=hamming(Nw,'periodic');     % omit sqrt if OV=4
    else
        W=sqrt(hamming(Nw,'periodic'));
    end
    W=W/sqrt(sum(W(1:Ns:Nw).^2));      % normalize window
    F_clean=rfft(enframe(s,W,Ns),Nw,2);      % do STFT: one row per time frame, +ve frequencies only
    F_noisy=rfft(enframe(x,W,Ns),Nw,2);
    
    % compute the true STFT noise spectrum (assumes additive noise distortion)
    noise = F_clean - F_noisy;
 
    % compute true SNR and threshold it to produce the ideal binary mask
    SNR = abs(F_clean).^2 ./ abs(noise).^2;     % IBM(t,f)=1 if Target(t,f)-Masker(t,f)>LC (dB)
    MASK = zeros( size(SNR) );
    MASK( SNR>10^(0.1*LC) ) = 1;        % set values to 1 if SNR of signal higher than threshold SNR
    
    % apply the ideal binary mask and create modified complex spectrum
    F_clean_IBM = abs(F_noisy) .* MASK .* exp(j*angle(F_noisy));
    
    y = overlapadd(irfft(F_clean_IBM,Nw,2),W,Ns);  % reconstitute the time waveform
    MASK = MASK';
    
% EOF
