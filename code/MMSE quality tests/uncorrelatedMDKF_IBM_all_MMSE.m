function Y_MDKF = uncorrelatedMDKF_IBM_all_MMSE(noisy, mmse_noisy, mmse_clean, fs, Tw, Ts, p, Tw_slow, Ts_slow, fs_slow, LC, u_present, var_present, u_absent, var_absent)

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

noisyfft = rfft(enframe(noisy,W,Ns),Nw,2);
noisy_fft = rfft(enframe(mmse_noisy,W,Ns),Nw,2);     % STFT of noisy signal
clean_fft = rfft(enframe(mmse_clean,W,Ns),Nw,2);

noisyfft_mag = abs(noisy_fft);
cleanfft_mag = abs(clean_fft);  % used to estimate LPCs in ideal case

f_noisy = enframe(mmse_noisy,W,Ns,'sp');
x = estnoiseg(f_noisy,Ns/fs);   % estimate the noise power spectrum; each row is 1 time frame

% after doing enframe, one row is one time frame (i.e. freq bins along rows)
% take transpose so that freq bins are along columns (easier to represent)
noisyfft_mag = noisyfft_mag';
noisy_fft = noisy_fft'; 
cleanfft_mag = cleanfft_mag';
noisyfft = noisyfft';
numBins = size(noisy_fft, 1);

% window and shift for each KF (in seconds)
[W_slow, Nw_slow, Ns_slow] = makeHammingWindow(fs_slow, Tw_slow, Ts_slow);

% mask used to scale mean and variance of observation
% enframe + rfft same dimensions as that applied in this function, so same numBins
[~, mask] = IdBM(mmse_noisy, mmse_clean, fs, Tw, Ts, LC);

filtered_matrix = zeros(size(noisyfft_mag));
for m = 1:numBins     % each frequency bin (rows of f) has its own KF
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
    
    varV = sum(x(:,m));
    P = varV*eye(p);      % error covariance matrix
    
    for j=1:size(noisyfft_mag,2)       % process each time sample individually within frequency bin
        % pick the LPC frame whose centre is closest to current index
        [~, closest] = min(abs(j*ones(size(centres)) - centres));   

        A(1,:) = -ar_coefs(closest, 2:p+1); 
        
        varW = energy_residual(closest)/length(cleanmag_frames(closest,:)); % variance of excitation, calculated from lpcauto

        % decorrelate observation from the rest of the state vector
        state = A*state;
        P = A*P*A' + varW*(d*d');
        
        % decompose covariance matrix into 4 blocks
        varObs = P(1,1);
        covarObs = P(2:p,1);    % first column of P, excluding top element
        % note that P(1,2:p) = covarObs'
        covarRest = P(2:p,2:p);     % remainder of P matrix
        
        % transform state vector: x -> z
        % z_n|n-1:  a priori state (predicted)
        T = eye(p);
        T(2:p,1) = -covarObs/varObs;
        z_predicted = T*state;
        z_Pmatrix = T*P*T';
        
        % multiply 3 gaussians together: mask (1 or 0), observation, prediction - 2 at a time
        % scale mean and variance of observation depending on whether mask gave 1 or 0 at this T-F unit
        if mask(m,j) == 1
            var1 = varV;
            var2 = var_present(m);
            var_y = 1/((1/var1 + 1/var2));
            u_y = var_y*(noisyfft_mag(m,j)/var1 + u_present(m)/var2);
        else
            var1 = varV;
            var2 = var_absent(m);
            var_y = 1/((1/var1 + 1/var2));
            u_y = var_y*(noisyfft_mag(m,j)/var1 + u_absent(m)/var2);
        end
        
        u_z = d'*z_predicted;
        var_z = varObs;
        
        % z_n|n:    a posteriori state (updated)
        var_zUpdated = (var_z*var_y)/(var_z + var_y);
        
        z_updated = z_predicted + var_z/(var_z+ var_y)*(u_y - u_z)*d;

        % get original state and covariance matrix
        state = T^(-1)*z_updated;
        
        z_Pmatrix(1,1) = var_zUpdated;
        P = T^(-1)*z_Pmatrix*((T^(-1))');
        
        filtered_matrix(m,j) = state(1);  % update estimated output
    end
    
    filtered_matrix(m,:) = filtered_matrix(m,:) .* exp(1i*angle(noisyfft(m,:)));
end

filtered_matrix = filtered_matrix';     % so that output time signal is col vector
Y_MDKF = overlapadd(irfft(filtered_matrix,Nw,2), W, Ns);
