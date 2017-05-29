function Y_MDKF = idealMDKF_linear(noisy, clean, fs, Tw, Ts, p)

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


% % initialisation for coloured noise
% d_v=zeros(q,1);
% d_v(1)=1;
% B=zeros(q);
% for b=1:q-1
%     A(b+1,b)=1;
% end


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


Tw_slow = 16e-3;     % window and shift for each KF
Ts_slow = 4e-3;
fs_slow = 1/Ts;
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

%         varV = sum(mean(x,1), 2);        % noise variance calculated from estnoiseg above
        
%         varV = x(m,j);
        
        varW = energy_residual(closest)/length(cleanmag_frames(closest,:)); % variance of excitation, calculated from lpcauto
    %     varW = var(cleanFrames(lpc_frame_index, :));

        state = A*state;
        P = A*P*A' + varW*(d*d');
        K = P*d*((varV + d'*P*d)^(-1));
        state = state + K*(noisyfft_mag(m,j) - d'*state);
        P = (eye(p) - K*d')*P;

        filtered_matrix(m,j) = state(1);  % update estimated output
    end
    
    filtered_matrix(m,:) = filtered_matrix(m,:) .* exp(1i*angle(noisy_fft(m,:)));
    
%     modulation_frames = enframe(noisyfft_mag(m,:), W_slow, Ns_slow);
%     angles = enframe(noisy_fft(m,:), W_slow, Ns_slow);
%     cleanmag_frames = enframe(cleanfft_mag(m,:), W_slow, Ns_slow);
    
%     state = modulation_frames(1,1:p)';  % initial state - doesn't really matter
    
    
    
%     for i = 1:size(modulation_frames, 1)        % number of frames
%         state = modulation_frames(i,1:p)';
%         
%         % LPCs and excitation variance constant within modulation frame
%         [ar, excitation_var] = lpcauto(cleanmag_frames(i,:),p);     % LPCs estimated from clean speech
% %         [ar, excitation_var] = lpcauto(modulation_frames(i,:),p);   % LPCs estimated from noisy speech
%         
%         A(1,:) = -ar(2:p+1);
%         
%         % estimate the noise power spectrum using current modulation frame
%         % comment out to use noise estimated from entire signal
%         x = estnoiseg((modulation_frames(i,:)),Ns_slow/fs_slow);
%         
%         % ------- subjective -------
%         % noise estimated from current frame results in less noisy output
%         % but either more musical noise or musical noise more obvious
%         varV = sum(mean(x,2),1);
%         
% %         varV = sum(x(:, 1 + (i-1)*Ns_slow),1);        % noise variance calculated using estnoiseg for each time frame
%         
%         varW = excitation_var/length(modulation_frames(i,:));
%         
%         P = varV*eye(p);              % initial error covariance matrix 
%         for j = 1:size(modulation_frames, 2)     % for each sample within modulation frame
%             state = A*state;
%             P = A*P*A' + varW*(d*d');
% 
%             K = P*d*((varV + d'*P*d)^(-1));
%             state = state + K*(modulation_frames(i,j) - d'*state);
%             P = (eye(p) - K*d')*P;
%             y(i,j) = state(1);    % update estimated output
%         end
%     end

%     f_filtered = y .* exp(1i*angle(angles));
%     filtered_matrix(m,:) = overlapadd(f_filtered, W_slow, Ns_slow);
end

filtered_matrix = filtered_matrix';     % so that output time signal is col vector
Y_MDKF = overlapadd(irfft(filtered_matrix,Nw,2), W, Ns);
