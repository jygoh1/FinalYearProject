function e = getExcitationVar(clean, W, Nw, Ns, p, W_slow, Ns_slow)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clean_fft = rfft(enframe(clean,W,Ns),Nw,2);

cleanfft_mag = abs(clean_fft);  % used to estimate LPCs in ideal case
cleanfft_mag = cleanfft_mag';
numBins = size(cleanfft_mag,1);

% window and shift for each KF (in seconds)

for m = 1:numBins     % each frequency bin (rows of f) has its own KF
    % assume |noisy| = |signal| + |noise| in modulation domain
    cleanmag_frames = enframe(cleanfft_mag(m,:), W_slow, Ns_slow);

    for i = 1:size(cleanmag_frames, 1)        % number of frames
        % LPCs and excitation variance constant within modulation frame
        [~, energy_residual(m,i)] = lpcauto(cleanmag_frames(i,:),p);     % LPCs estimated from clean speech
        energy_residual(m,i) = energy_residual(m,i)/length(cleanmag_frames(i,:));
    end
end

e = mean(mean(energy_residual));    % average residual excitation energy over all frames

end

