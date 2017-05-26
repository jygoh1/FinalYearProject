function [ y, MASK, SNR ] = idbm( x, s, fs, Tw, Ts, LC )
% IDBM Apply spectral ideal binary mask to noisy speech.
%
%   [Y,MASK,SNR]=IDBM(X,S,FS,TW,TS,LC) performs FFT based short-time
%   spectral analysis-modification-synthesis, with frame duration 
%   of TW ms and frame shift of TS ms. Spectral modification 
%   performed is the application of an ideal binary mask to noise
%   corrupted speech signal. Additive noise distortion is assumed.
%   The noisy speech signal is given in vector X and the 
%   corresponding clean speech signal is given in vector S. 
%   Each signal is sampled at FS Hz. The ideal binary mask is 
%   computed from an oracle (true) signal-to-noise ratio (SNR) by
%   thresholding with local SNR criterion specified in LC. 
%   The synthesized enhanced speech is returned in Y, the ideal 
%   binary mask is returned in MASK, while the true instantaneous
%   spectral signal-to-noise ratio is returned in SNR.
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
%           SNR is the true (oracle) instantaneous spectral SNR (dB).
%
%   Example
%           % read clean speech samples from wav file
%           [ clean, fs, nbits ] = wavread( 'sp10.wav' ); 
%       
%           % read noise samples from wav file
%           [ noise, fs, nbits ] = wavread( 'babble.wav' ); 
%       
%           % mix target speech and the masker noise at 0 dB SNR
%           noisy = addnoise( clean, noise, 0 ); 
%       
%           Tw = 32;    % analysis frame duration (ms) 
%           Ts = Tw/8;  % analysis frame shift (ms)
%           LC = -5;    % local SNR criterion (dB)
%       
%           % apply ideal binary mask to noise corrupted speech signal
%           enhanced = idbm( noisy, clean, fs, 32, 32/8, -5 );
%       
%           % plot all signal 
%           figure('Position', [ 25 100 800 600 ], 'PaperPositionMode', 'auto', 'color', 'w' );
%           time = [ 0:length(noisy)-1 ]/fs;
%           
%           subplot( 3,1,1 );
%           plot( time, clean, 'k-' );
%           xlim( [ min(time) max(time) ] );
%           title( sprintf('Clean speech') ); 
%           xlabel( sprintf('Time (s)') ); 
%           ylabel( sprintf('Amplitude') ); 
%       
%           subplot( 3,1,2 );
%           plot( time, noisy, 'k-' );
%           xlim( [ min(time) max(time) ] );
%           title( sprintf('Noisy speech') ); 
%           xlabel( sprintf('Time (s)') ); 
%           ylabel( sprintf('Amplitude') ); 
%       
%           subplot( 3,1,3 );
%           plot( time, enhanced, 'k-' );
%           xlim( [ min(time) max(time) ] );
%           title( sprintf('Ideal binary mask enhanced speech') ); 
%           xlabel( sprintf('Time (s)') ); 
%           ylabel( sprintf('Amplitude') ); 
%       
%           print( '-dpng', 'idbm.png'  ); 
%
%   See also TEST_IDBM.

%   Author: Kamil Wojcicki, UTD, October 2011

    
    % check for correct number of input arguments
    if( nargin~=6 ), error( sprintf('Not enough input arguments. Type "help %s" for usage help.', mfilename) ); end;

    L = length( x );                        % length of the speech signal
    Nw = round( fs*Tw*0.001 );              % frame duration (in samples)
    Ns = round( fs*Ts*0.001 );              % frame shift (in samples)
    nfft = 2^nextpow2( 2*Nw );              % FFT analysis length

    % divide noisy and clean speech signals into frames
    [ frames.x, indexes ] = vec2frames( x, Nw, Ns, 'rows', @hanning, true );
    [ frames.s, indexes ] = vec2frames( s, Nw, Ns, 'rows', @hanning, true );

    % perform short-time Fourier transform (STFT) analyses
    spectrum.X = fft( frames.x, nfft, 2 );           
    spectrum.S = fft( frames.s, nfft, 2 );

    % compute the true STFT noise spectrum (assumes additive noise distortion)
    spectrum.D = spectrum.X - spectrum.S;
 
    % compute true SNR and threshold it to produce the ideal binary mask
    SNR = abs(spectrum.S).^2 ./ abs(spectrum.D).^2;
    MASK = zeros( size(SNR) );
    MASK( SNR>10^(0.1*LC) ) = 1;

    % apply the ideal binary mask and create modified complex spectrum
    spectrum.Y = abs(spectrum.X) .* MASK .* exp(j*angle(spectrum.X));

    % alternatively, you could perform a sanity check and reconstruct the clean speech...
    % spectrum.Y = abs(spectrum.S) .* exp(j*angle(spectrum.S));

    % apply inverse STFT 
    frames.y = real( ifft(spectrum.Y,nfft,2) ); 

    % discard FFT padding from frames
    frames.y = frames.y(:, 1:Nw);

    % perform overlap-and-add synthesis
    y = frames2vec( frames.y, indexes, 'rows', @hanning, 'G&L' ); 

    % truncate extra padding (back to original signal length)
    y = y(1:L);


% EOF
