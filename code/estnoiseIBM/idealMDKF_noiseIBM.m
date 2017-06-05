function Y_MDKF = idealMDKF_noiseIBM(noisy, clean, fs, Tw, Ts, p, Tw_slow, Ts_slow, fs_slow, LC, mask_th)

%   noisy: noisy input speech
%   clean: clean input speech
%   fs: sampling frequency
%   Tw: frame duration in s
%   Ts: frame shift in s (overlap)
%   p: MDKF order
%   q: (for coloured noise case)
%
%   MDKF: Apply modulation-domain Kalman filtering for speech enhancement
%   Author: Jia Ying Goh, Imperial College London, January 2017
%   Credits:
%       Mike Brookes, Imperial College London - VOICEBOX: A speech processing toolbox for MATLAB

d=zeros(p,1);
d(1)=1;

% state transition matrix
A=zeros(p);
for a=1:p-1
    A(a+1,a)=1;
end

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

noisy_fft = rfft(enframe(noisy,W,Ns),Nw,2);     % STFT of noisy signal
clean_fft = rfft(enframe(clean,W,Ns),Nw,2);

noisyfft_mag = abs(noisy_fft);
cleanfft_mag = abs(clean_fft);  % used to estimate LPCs in ideal case

% mask to modify noise estimation
[~, mask] = IdBM(noisy, clean, fs, Tw, Ts, LC);

f_noisy = enframe(noisy,W,Ns,'sp');
x = estnoiseIBM(f_noisy,Ns/fs,mask',mask_th);   % estimate the noise power spectrum; each row is 1 time frame
noisepow = 10*log10(sum(mean(x,1),2))

% after doing enframe, one row is one time frame (i.e. freq bins along rows)
% take transpose so that freq bins are along columns (easier to represent)
noisyfft_mag = noisyfft_mag';
noisy_fft = noisy_fft'; 
cleanfft_mag = cleanfft_mag';

% window and shift for each KF (in seconds)
[W_slow, Nw_slow, Ns_slow] = makeHammingWindow(fs_slow, Tw_slow, Ts_slow);

filtered_matrix = zeros(size(noisyfft_mag));

for m = 1:size(noisyfft_mag,1)     % each frequency bin (rows of f) has its own KF
    % assume |noisy| = |signal| + |noise| in modulation domain
    cleanmag_frames = enframe(cleanfft_mag(m,:), W_slow, Ns_slow);
    
    numFrames = size(cleanmag_frames, 1);
    framelen = Nw_slow;      % length of each frame

    centres = zeros(numFrames, 1);
    centres(1) = round(framelen/2);
    for i = 2:numFrames
        centres(i) = round(i*Ns_slow + framelen/2);
    end

    for i = 1:size(cleanmag_frames, 1)        % number of frames
        % LPCs and excitation variance constant within modulation frame
        [ar_coefs(i,:), energy_residual(i)] = lpcauto(cleanmag_frames(i,:),p);     % LPCs estimated from clean speech
    end
    
    state = noisyfft_mag(m,1:p)';  % initial state - not that important
    
%     varV = sum(mean(x,1),2);
    varV = sum(x(:,m));
    P = varV*eye(p);      % error covariance matrix
    
    for j=1:size(noisyfft_mag,2)       % process each time sample individually within frequency bin
        % pick the LPC frame whose centre is closest to current index
        [~, closest] = min(abs(j*ones(size(centres)) - centres));   

        A(1,:) = -ar_coefs(closest, 2:p+1);
        
        varW = energy_residual(closest)/length(cleanmag_frames(closest,:)); % variance of excitation, calculated from lpcauto

        state = A*state;
        P = A*P*A' + varW*(d*d');
        K = P*d*((varV + d'*P*d)^(-1));
        state = state + K*(noisyfft_mag(m,j) - d'*state);
        P = (eye(p) - K*d')*P;

        filtered_matrix(m,j) = state(1);  % update estimated output
    end
    
    filtered_matrix(m,:) = filtered_matrix(m,:) .* exp(1i*angle(noisy_fft(m,:)));
end

filtered_matrix = filtered_matrix';     % so that output time signal is col vector
Y_MDKF = overlapadd(irfft(filtered_matrix,Nw,2), W, Ns);
