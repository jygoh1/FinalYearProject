function Y_MDKF = ExptIbmKF(noisy, clean, fs, Tw, Ts, p, LC)

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



% time-domain signal -> STFT, each modulation signal is time signal for one freq bin over time
% each noisy modulation signal |Y(n,k)| windowed into short modulation frames
% LPCs and excitation var estimated for each frame
% perform the KF for each freq bin (in each column, using enframe)
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

f_noisy = enframe(noisy,W,Ns,'sp');
x = estnoiseg(f_noisy,Ns/fs);   % estimate the noise power spectrum; each row is 1 time frame
noisepow = 10*log10(sum(mean(x,1),2));
noisepow

% after doing enframe, one row is one time frame (i.e. freq bins along rows)
% take transpose so that freq bins are along columns (easier to represent)
noisyfft_mag = noisyfft_mag';
noisy_fft = noisy_fft'; 
cleanfft_mag = cleanfft_mag';


numBins = size(noisy_fft, 1);

% get binary mask - each row is one freq bin
[~, mask] = IdBM(noisy, clean, fs, Tw, Ts, LC);
present = mask .* noisyfft_mag;
absent = noisyfft_mag - present;

u_present = zeros(numBins, 1);
var_present = zeros(numBins, 1);
u_absent = zeros(numBins, 1);
var_absent = zeros(numBins, 1);
for i = 1:numBins
    present_elements = (present(i,(present(i,:) > 0)));
    u_present(i) = mean(present_elements);    % mean of each frequency bin
    var_present(i) = var(present_elements);
    
    absent_elements = (absent(i,(absent(i,:) > 0)));
    u_absent(i) = mean(absent_elements);    % mean of each frequency bin
    var_absent(i) = var(absent_elements);
end


for m = 1:numBins           % each frequency bin (rows of f) has its own KF
    Tw_slow = 16e-3;        % window and shift for each KF
    Ts_slow = 4e-3;
    fs_slow = 1/Ts;
    [W_slow, Nw_slow, Ns_slow] = makeHammingWindow(fs_slow, Tw_slow, Ts_slow);
    
    modulation_frames = enframe(noisyfft_mag(m,:), W_slow, Ns_slow);
    angles = enframe(noisy_fft(m,:), W_slow, Ns_slow);      % noisy phase used as phase of output
    cleanmag_frames = enframe(cleanfft_mag(m,:), W_slow, Ns_slow);
    
    % assume |noisy| = |signal| + |noise| in modulation domain
    
    state = modulation_frames(1, 1:p)';  % initial state - doesn't really matter
    
    y = zeros(size(modulation_frames));
    
    for i = 1:size(modulation_frames, 1)        % number of frames
        % LPCs and excitation variance constant within modulation frame
        [ar, excitation_var] = lpcauto(cleanmag_frames(i,:), p);     % LPCs estimated from clean speech
%         [ar, excitation_var] = lpcauto(modulation_frames(i,:), p);   % LPCs estimated from noisy speech
        
        A(1,:) = -ar(2:p+1);
        
        % estimate the noise power spectrum using current modulation frame
        % comment out to use noise estimated from entire signal
        x = estnoiseg((modulation_frames(i,:)), Ns_slow/fs_slow);
        
        % ------- subjective -------
        % noise estimated from current frame results in less noisy output
        % but either more musical noise or musical noise more obvious
        varV = sum(mean(x,2),1);
        
        varW = excitation_var/length(modulation_frames(i,:));
        
        P = varV*eye(p);              % initial error covariance matrix 
        for j = 1:size(modulation_frames, 2)     % for each sample within modulation frame
            state = A*state;
            P = A*P*A' + varW*(d*d');

            K = P*d*((varV + d'*P*d)^(-1));
            state = state + K*(modulation_frames(i,j) - d'*state);
            P = (eye(p) - K*d')*P;
            y(i,j) = state(1);    % update estimated output
        end
    end

    f_filtered = y .* exp(1i*angle(angles));
    filtered_matrix(m,:) = overlapadd(f_filtered, W_slow, Ns_slow);
end

filtered_matrix = filtered_matrix';     % so that output time signal is col vector
Y_MDKF = overlapadd(irfft(filtered_matrix,Nw,2), W, Ns);
