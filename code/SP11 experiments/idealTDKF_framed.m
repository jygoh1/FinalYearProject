function Y_TDKF = idealTDKF_framed(noisy, clean, fs, Tw, Ts, p)

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
% ar = lpcauto(clean,p);
% A(1,:) = -ar(2:p+1);

state = noisy(1:p);       % initial state - doesn't seem to have any (significant) effect on output
% varV = var(noisy);
% P = varV*eye(p);      % error covariance matrix

[W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts);

noisyFrames = enframe(noisy,W,Ns);                          % do processing for each time frame

cleanFrames = enframe(clean,W,Ns);

% ar = lpcauto(cleanFrames(1,:),p);
% A(1,:) = -ar(2:p+1);

%%

% need 'sp' to scale window so that adding power in all +ve freq DFT bins gives the total power in the original signal
f_noisy = enframe(noisy, W, Ns, 'sp');
x = estnoiseg(f_noisy, Ns/fs); % estimate the noise power spectrum
noisepow = 10*log10(sum(mean(x,1), 2));
noisepow
% figure;
% plot(noiselev,pow,noiselev,noiselev,':k');
% xlabel('True noise power (dB) [speech = 0 dB]');
% ylabel('Estimated Noise Power (dB)');
% legend('white','location','northwest');

varV = sum(mean(x,1),2);
P = varV*eye(p);      % error covariance matrix


%%

for i = 1:size(cleanFrames, 1)
    [ar_coefs(i,:), energy_residual(i)] = lpcauto(cleanFrames(i,:), p);   % estimate LPC coefficients from clean or noisy speech (each frame)
end

numFrames = size(noisyFrames, 1);
framelen = size(noisyFrames, 2);      % length of each frame
cutoff = round(framelen/4);


% Y_TDKF = noisy;
% for i=1:length(noisy)       % process each time sample individually
%     % pick LPC frame for each sample as close to centre of frame as possible
%     % first 1/4 frame belongs to A_1, last 1/4 frame belongs to A_k
%     % all others: pick frame index by (i - 1/4 frame) / framelength
%     lpc_frame_index = max(round((i-cutoff)/framelen) + 1, 1);
%     
%     lpc_frame_index = max(round(i/framelen) + 1, 1);
%     
%     lpc_frame_index = min(lpc_frame_index, numFrames);
%     A(1,:) = -ar_coefs(lpc_frame_index, 2:p+1);  
%     
%     varV = sum(mean(x,1), 2);        % noise variance calculated from estnoiseg above
%     varW = energy_residual(lpc_frame_index); % variance of excitation, calculated from lpcauto
% %     varW = var(cleanFrames(lpc_frame_index, :));
%     
%     state = A*state;
%     P = A*P*A' + varW*(d*d');
% %             P = A*P*A' + varW*eye(p);
% 
%     K = P*d*((varV + d'*P*d)^(-1));
%     state = state + K*(noisy(i) - d'*state);
%     P = (eye(p) - K*d')*P;
%     Y_TDKF(i) = state(1);  % update estimated output
% end


Y_TDKF = noisyFrames;
for i=1:numFrames
    A(1,:) = -ar_coefs(i, 2:p+1);
    
%     varV=var(noisyFrames(i,:));  
    varV = sum(mean(x,1),2);        % noise variance calculated from estnoiseg above
    
%     varW = var(cleanFrames(i,:));            % this seems to work better than using energy_residual(i)?
    varW = energy_residual(i)/length(cleanFrames(i,:));
    
    for j=1:size(noisyFrames,2)   % number of data points in each frame; all other variables are updated per sample
        state = A*state;
        P = A*P*A' + varW*(d*d');
%         P = A*P*A' + varW*eye(p);

        K = P*d*((varV + d'*P*d)^(-1));
        state = state + K*(noisyFrames(i,j) - d'*state);
        P = (eye(p) - K*d')*P;
        Y_TDKF(i,j) = state(1);  % update estimated output
    end
end

Y_TDKF = overlapadd(Y_TDKF,W,Ns);


% EOF
