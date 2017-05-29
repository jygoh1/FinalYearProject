function Y_TDKF = idealTDKF_linear(noisy, clean, fs, Tw, Ts, p)

%   noisy: noisy input speech
%   clean: clean input speech
%   fs: sampling frequency
%   Tw: frame duration in s
%   Ts: frame shift in s (overlap)
%   p: TDKF order
%
%   TDKF: Apply time-domain Kalman filtering for speech enhancement
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

state = noisy(1:p);       % initial state - doesn't seem to have any (significant) effect on output

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

noisyFrames = enframe(noisy,W,Ns);                          % do processing for each time frame

cleanFrames = enframe(clean,W,Ns);

%%

% need 'sp' to scale window so that adding power in all +ve freq DFT bins gives the total power in the original signal
f_noisy = enframe(noisy, W, Ns, 'sp');
x = estnoiseg(f_noisy, Ns/fs); % estimate the noise power spectrum
noisepow = 10*log10(sum(mean(x,1), 2));
noisepow

varV = sum(mean(x,1),2);
P = varV*eye(p);      % error covariance matrix

%%

for i = 1:size(cleanFrames, 1)
    [ar_coefs(i,:), energy_residual(i)] = lpcauto(cleanFrames(i,:), p);   % estimate LPC coefficients from clean or noisy speech (each frame)
end

numFrames = size(noisyFrames, 1);
framelen = Nw;      % length of each frame

centres = zeros(numFrames, 1);
centres(1) = round(framelen/2);
for i = 2:numFrames
    centres(i) = round(i*Ns + framelen/2);
end

Y_TDKF = zeros(size(noisy));
for i=1:length(noisy)       % process each time sample individually
    % pick the LPC frame whose centre is closest to current index
    [~, closest] = min(abs(i*ones(size(centres)) - centres));   
    
    A(1,:) = -ar_coefs(closest, 2:p+1);  
    
    varV = sum(mean(x,1), 2);        % noise variance calculated from estnoiseg above

    varW = energy_residual(closest)/length(cleanFrames(closest,:)); % variance of excitation, calculated from lpcauto
%     varW = var(cleanFrames(lpc_frame_index, :));
    
    state = A*state;
    P = A*P*A' + varW*(d*d');
    K = P*d*((varV + d'*P*d)^(-1));
    state = state + K*(noisy(i) - d'*state);
    P = (eye(p) - K*d')*P;
    
    Y_TDKF(i) = state(1);  % update estimated output
end

