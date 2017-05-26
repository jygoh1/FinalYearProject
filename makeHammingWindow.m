function [W, Nw, Ns] = makeHammingWindow(fs, Tw, Ts)

% fs: sampling frequency
% Tw: frame duration in s
% Ts: frame shift in s
% 
% W: output window
% Nw: frame duration in number of samples
% Ns: frame shift in number of samples

Nw = round(fs*Tw);                      % frame duration (in samples)
Ns = round(fs*Ts);                      % frame shift (in samples)
if Nw/Ns==4
    W=hamming(Nw,'periodic');                   % omit sqrt if OV=4
else
    W=sqrt(hamming(Nw,'periodic'));
end
W = W/sqrt(sum(W(1:Ns:Nw).^2));                   % normalize window

end

