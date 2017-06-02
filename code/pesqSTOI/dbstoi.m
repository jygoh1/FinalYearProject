function val = dbstoi(xl, xr, yl, yr, fs)
% -----------------------------------------------------------------------------%
%                               Inputs                                         %
% -----------------------------------------------------------------------------%
% xl : Clean speech signal from left ear.
% xr : Clean speech signal from right ear.
% yl : Noisy/processed signal from left ear.
% yr : Noisy/processed signal from right ear.
% fs : Common sampling frequency of the input signals [Hz].
%
% -----------------------------------------------------------------------------%
%                               Description                                    %
% -----------------------------------------------------------------------------%
% This is an implementation of the Deterministic Binaural Short-Time
% Objective Intelligibility (DBSTOI) measure as described in:
%
%   A. H. Andersen, J. M. de Haan, Z.-H. Tan, and J. Jensen, “Predicting the
%   Intelligibility of Noisy and Non-Linearly Processed Binaural Speech,”
%   Transactions on Audio Speech and Language Processing, <...>. 
% 
% The code is commented with references to equation numbers whenever
% possible. We simply use "STOI" to refer to the modified version of the STOI 
% measure which is part of the DBSTOI measure. 
%
% - Asger Heidemann Andersen, 05/07-2016
%
% -----------------------------------------------------------------------------%
%                                   Licence                                    %
% -----------------------------------------------------------------------------%
% License Agreement for the DBSTOI software (the ”Software”)
% 
% By using the Software you agree to the below terms and conditions. If you do 
% not agree with such terms and conditions, do not use the software.
% 
% All title, copyrights and pending patents in and to the Software are owned by 
% Oticon A/S and/or Aalborg University (the “Owners”). 
% 
% The Software may only be used for non-commercial purposes such as research 
% undertaken by public or private institutions and/or universities. Hence, the 
% Software may not be used in the context of product development.
% 
% DISCLAIMER OF WARRANTIES: THE SOFTWARE IS BEING PROVIDED TO YOU "AS IS" 
% WITHOUT WARRANTY, UPGRADES OR SUPPORT OF ANY KIND. THE OWNERS EXPRESSLY 
% DISCLAIM ALL WARRANTIE WITH REGARD TO THE SOFTWARE, EXPRESS OR IMPLIED, 
% INCLUDING, WITHOUT LIMITATION, ANY IMPLIED WARRANTIES OF FITNESS FOR A 
% PARTICULAR PURPOSE, QUALITY, ACCURACY OR NON-INFRINGEMENT OF THIRD PARTY 
% RIGHTS. 
% 
% LIMITATION OF LIABILITY: IN NO EVENT THE OWNERS WILL BE LIABLE TO YOU FOR ANY 
% INCIDENTIAL, SPECIAL, INDIRECT, CONSEQUENTIAL, OR PUNITIVE DAMAGES ARISING 
% OUT OF OR IN CONNECTION WITH YOUR USE OF THE SOFTWARE, WHETHER UNDER A THEORY 
% OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCTS LIABILITY, OR 
% OTHERWISE, EVEN IF THE OWNERS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH 
% DAMAGES. 


%% ----------------------------------------------------------------------------%
%                       Parameters and Constants                               %
% -----------------------------------------------------------------------------%
% Parameters for STOI
fs_int      = 10000;    % Sample rate of proposed intelligibility measure [Hz].
N_frame    	= 256;      % Window support [samples].
K           = 512;      % FFT size [samples].
J           = 15;       % Number of 1/3 octave bands [-].
mn          = 150;      % Center frequency of first 1/3 octave band [Hz].
N           = 30;       % Number of frames for intermediate intelligibility 
                        % measure (length analysis window) [-].
dyn_range   = 40;       % Speech dynamic range [dB].

% Parameters for EC-parameters search
tau_min     = -0.001;   % Minumum interaural delay compensation [s]. 
tau_max     = 0.001;    % Maximum interaural delay compensation [s].
ntaus       = 100;      % Number of tau values to try out [-]. 
gamma_min   = -20;      % Minumum interaural level compensation [dB].
gamma_max   = 20;       % Maximum interaural level compensation [dB].
ngammas     = 40;       % Number of gamma values to try out [-].

% Constants for jitter in EC-stage
sig_del_0   = 65e-6;    % ITD compensation standard deviation [s].
sig_eps_0   = 1.5;      % ILD compensation standard deviation [-].
alpha_0_db  = 13;       % Constant for level shift deviation [dB].
tau_0       = 1.6e-3;   % Constant for time shift deviation [s].
p           = 1.6;      % Constant for level shift deviation [-].

%% ----------------------------------------------------------------------------%
%                           Signal Preparation                                 %
% -----------------------------------------------------------------------------%
% Ensure that inputs are column vectors.
xl = xl(:);
xr = xr(:);
yl = yl(:);
yr = yr(:);

% Resample signals to internal sampling rate.
xl = resample(xl, fs_int, fs);
xr = resample(xr, fs_int, fs);
yl = resample(yl, fs_int, fs);
yr = resample(yr, fs_int, fs);

% Remove silent frames.
[xl, xr, yl, yr] = ...
    remove_silent_frames(xl, xr, yl, yr, dyn_range, N_frame, N_frame/2);
                                              
%% ----------------------------------------------------------------------------%
%                               STDFT and filtering                            %
% -----------------------------------------------------------------------------%
% Get 1/3 octave band matrix.
[H, cf, fids] = thirdoct(fs, K, J, mn);  	
cf = 2*pi*cf; 	

% Apply short time DFT to signals.
xl_hat = stdft(xl, N_frame, N_frame/2, K);
xr_hat = stdft(xr, N_frame, N_frame/2, K);
yl_hat = stdft(yl, N_frame, N_frame/2, K);
yr_hat = stdft(yr, N_frame, N_frame/2, K);

% Take single sided spectrum of signals.
xl_hat = xl_hat(:, 1:(K/2+1)).';
xr_hat = xr_hat(:, 1:(K/2+1)).';
yl_hat = yl_hat(:, 1:(K/2+1)).';
yr_hat = yr_hat(:, 1:(K/2+1)).';

%% ----------------------------------------------------------------------------%
%                   Compute intermediate correlation by (iii)                  %
% -----------------------------------------------------------------------------%
% Grid for intermediate correlation.
d = zeros(J, length(N:size(xl_hat, 2)));

% Interaural compensation time and delay values.
taus = linspace(tau_min, tau_max, ntaus);
gammas = linspace(gamma_min, gamma_max, ngammas);
sigma_epsilon = sqrt(2)*sig_eps_0*(1+(abs(gammas)./alpha_0_db).^p)/20;
gammas = gammas/20;
sigma_delta = sqrt(2) * sig_del_0 * (1 + (abs(taus)/tau_0));       

% Evaluate (11) for all third octave bands, all frames, and all searched
% gamma/tau combinations. This involves evaluating (12) as well as
% corresponding equations for the denominator of (11).
for i=1:J % Loop 1/3 octave bands.
    for j=1:size(d,2) % Loop frames.
        seg_xl = xl_hat(fids(i,1):fids(i,2), j:(j+N-1));
        seg_xr = xr_hat(fids(i,1):fids(i,2), j:(j+N-1));
        seg_yl = yl_hat(fids(i,1):fids(i,2), j:(j+N-1));
        seg_yr = yr_hat(fids(i,1):fids(i,2), j:(j+N-1));

        Lx = sum(seg_xl'.'.*seg_xl);
        Lx = Lx - mean(Lx);
        Rx = sum(seg_xr'.'.*seg_xr);
        Rx = Rx - mean(Rx);
        rhox = sum(seg_xr'.'.*seg_xl);
        rhox = rhox - mean(rhox);
        Ly = sum(seg_yl'.'.*seg_yl);
        Ly = Ly - mean(Ly);
        Ry = sum(seg_yr'.'.*seg_yr);
        Ry = Ry - mean(Ry);
        rhoy = sum(seg_yr'.'.*seg_yl);
        rhoy = rhoy - mean(rhoy); 

        % Evaluate parts of intermediate correlation.
        exy = ...
            ones(ntaus,1) * ((10.^(2*gammas)*sum(Lx.*Ly) + 10.^(-2*gammas)*sum(Rx.*Ry)) .* exp(2*log(10)^2 .* sigma_epsilon.^2)) + ...
            sum(Lx.*Ry) + sum(Rx.*Ly) - ...
            2*((Lx*real(rhoy.'*exp(-1i*cf(i)*taus))+Ly*real(rhox.'*exp(-1i*cf(i)*taus)))'*10.^(gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) - ...
            2*((Rx*real(rhoy.'*exp(-1i*cf(i)*taus))+Ry*real(rhox.'*exp(-1i*cf(i)*taus)))'*10.^(-gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) + ...
            2*(real(rhox*rhoy') + exp(-2*cf(i)^2*sigma_delta.^2).*real(rhox*rhoy.'*exp(-1i*2*cf(i)*taus)))'*ones(1,ngammas);     
        exx = ...
            ones(ntaus,1) * ((10.^(2*gammas)*sum(Lx.*Lx) + 10.^(-2*gammas)*sum(Rx.*Rx)) .* exp(2*log(10)^2 .* sigma_epsilon.^2)) + ...
            2*sum(Lx.*Rx) - ...
            2*((Lx*real(rhox.'*exp(-1i*cf(i)*taus))+Lx*real(rhox.'*exp(-1i*cf(i)*taus)))'*10.^(gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) - ...
            2*((Rx*real(rhox.'*exp(-1i*cf(i)*taus))+Rx*real(rhox.'*exp(-1i*cf(i)*taus)))'*10.^(-gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) + ...
            2*(real(rhox*rhox') + exp(-2*cf(i)^2*sigma_delta.^2).*real(rhox*rhox.'*exp(-1i*2*cf(i)*taus)))'*ones(1,ngammas);
        eyy = ...
            ones(ntaus,1) * ((10.^(2*gammas)*sum(Ly.*Ly) + 10.^(-2*gammas)*sum(Ry.*Ry)) .* exp(2*log(10)^2 .* sigma_epsilon.^2)) + ...
            2*sum(Ly.*Ry) - ...
            2*((Ly*real(rhoy.'*exp(-1i*cf(i)*taus))+Ly*real(rhoy.'*exp(-1i*cf(i)*taus)))'*10.^(gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) - ...
            2*((Ry*real(rhoy.'*exp(-1i*cf(i)*taus))+Ry*real(rhoy.'*exp(-1i*cf(i)*taus)))'*10.^(-gammas)).*exp(0.5*(ones(ntaus,1)*(log(10))^2*sigma_epsilon.^2 - cf(i)^2*sigma_delta'.^2*ones(1,ngammas))) + ...
            2*(real(rhoy*rhoy') + exp(-2*cf(i)^2*sigma_delta.^2).*real(rhoy*rhoy.'*exp(-1i*2*cf(i)*taus)))'*ones(1,ngammas);  

        % Evaluate expected intermediate correlation.
        d(i,j) = max(max(exy./sqrt(exx.*eyy)));

    end
end

%% ----------------------------------------------------------------------------%
%                           Compute better ear STOI                            %
% -----------------------------------------------------------------------------%
% Arrays for the 1/3 octave envelope.
Xl = zeros(J, size(xl_hat, 2));
Xr = zeros(J, size(xl_hat, 2));
Yl = zeros(J, size(xl_hat, 2));
Yr = zeros(J, size(xl_hat, 2));

% Apply 1/3 octave bands.
for k=1:size(xl_hat,2) % Loop frames.
    Xl(:, k)	= H*abs(xl_hat(:, k)).^2;
    Xr(:, k)	= H*abs(xr_hat(:, k)).^2;
    Yl(:, k)	= H*abs(yl_hat(:, k)).^2;
    Yr(:, k)	= H*abs(yr_hat(:, k)).^2;
end

% Arrays for better-ear correlations.
dl_interm  	= zeros(J, length(N:size(xl_hat, 2)));
dr_interm  	= zeros(J, length(N:size(xl_hat, 2)));

% Compute intermediate better-ear correlations.
for m = N:size(xl_hat,2) % Loop frames.
    Xl_seg  	= Xl(:, (m-N+1):m);
    Xr_seg  	= Xr(:, (m-N+1):m);
    Yl_seg  	= Yl(:, (m-N+1):m);
    Yr_seg  	= Yr(:, (m-N+1):m);
    for n = 1:J % Loop 1/3 ocatve bands.
        xln = Xl_seg(n,:).'-sum(Xl_seg(n,:))/N;
        xrn = Xr_seg(n,:).'-sum(Xr_seg(n,:))/N; 
        yln = Yl_seg(n,:).'-sum(Yl_seg(n,:))/N;
        yrn = Yr_seg(n,:).'-sum(Yr_seg(n,:))/N;
        dl_interm(n,m-N+1) = sum(xln.*yln)/(norm(xln)*norm(yln));
        dr_interm(n,m-N+1) = sum(xrn.*yrn)/(norm(xrn)*norm(yrn));
    end
end

% Get the better ear intermediate coefficients.
dl_interm(or(isnan(dl_interm),isinf(dl_interm))) = 0;
dr_interm(or(isnan(dr_interm),isinf(dr_interm))) = 0;
dbe_interm = max(dl_interm, dr_interm);

%% ----------------------------------------------------------------------------%
%                               Compute STOI measure                           %
% -----------------------------------------------------------------------------%
d = max(d, dbe_interm);
val = mean(d(:));

%% ----------------------------------------------------------------------------%
%                               Auxiliary functions                            %
% -----------------------------------------------------------------------------%
function [xl_sil, xr_sil, yl_sil, yr_sil] = ...
                remove_silent_frames(xl, xr, yl, yr, range, N, K)
% Remove region in signals where the target speaker is silent.
% Inputs:
%   xl:         Clean speech signal from left ear.
%   xr:         Clean speech signal from right ear.
%   nl:         Noisy/processed signal from left ear..
%   nr:         Noisy/processed signal from right ear..
%   range:      Dynamic range of speech [dB].
%   N:          Length of frames [samples].
%   K:          Overlap of frames [samples].
% Outputs:
%   xl_sil:     Clean signal, left ear, with only active speech.
%   xr_sil:     Clean signal, right ear, with only active speech.
%   yl_sil:     Noisy/processed signal, left ear, with only active speech.
%   yr_sil:     Noisy/processed signal, right ear, with only active speech.

xl       = xl(:); 
xr       = xr(:);
yl       = yl(:);
yr       = yr(:);

frames  = 1:K:(length(xl)-N);
w       = hanning(N);
msk_l     = zeros(size(frames));
msk_r     = zeros(size(frames));

for j = 1:length(frames)
    jj      = frames(j):(frames(j)+N-1);
    msk_l(j) 	= 20*log10(norm(xl(jj).*w)./sqrt(N));
    msk_r(j) 	= 20*log10(norm(xr(jj).*w)./sqrt(N));
end

msk_l   = (msk_l-max([msk_l msk_r])+range)>0;
msk_r   = (msk_r-max([msk_l msk_r])+range)>0;
msk     = or(msk_l, msk_r);
count   = 1;

xl_sil   = zeros(size(xl));
xr_sil   = zeros(size(xr));
yl_sil   = zeros(size(yl));
yr_sil   = zeros(size(yr));

for j = 1:length(frames)
    if msk(j)
        jj_i            = frames(j):(frames(j)+N-1);
        jj_o            = frames(count):(frames(count)+N-1);
        xl_sil(jj_o)    = xl_sil(jj_o) + xl(jj_i).*w;
        xr_sil(jj_o)    = xr_sil(jj_o) + xr(jj_i).*w;
        yl_sil(jj_o)  	= yl_sil(jj_o) + yl(jj_i).*w;
        yr_sil(jj_o)  	= yr_sil(jj_o) + yr(jj_i).*w;
        count           = count+1;
    end
end

xl_sil = xl_sil(1:jj_o(end));
xr_sil = xr_sil(1:jj_o(end));
yl_sil = yl_sil(1:jj_o(end));
yr_sil = yr_sil(1:jj_o(end));

end


function  [A, cf, fids] = thirdoct(fs, n_fft, num_bands, mn)
% Generate 1/3 octave matrix.
% Inputs:
%   fs:         Samplerate [Hz].
%  	n_fft:      FFT size [samples].
%  	num_bands:  Number of bands [-].
%  	mn:         Center frequency of first 1/3 octave band [Hz].
% Outputs:
% 	A:          Octave band matrix [-].
% 	cf:         Center frequencies [Hz].
%  	fids:       Indices of frequency edges [-].

f               = linspace(0, fs, n_fft+1);
f               = f(1:(n_fft/2+1));
k               = 0:(num_bands-1); 
cf              = 2.^(k/3)*mn;
fl              = sqrt((2.^(k/3)*mn).*2.^((k-1)/3)*mn);
fr              = sqrt((2.^(k/3)*mn).*2.^((k+1)/3)*mn);
A               = zeros(num_bands, length(f));
fids            = zeros(num_bands,2); 

for i = 1:(length(cf))
    [a b]                   = min((f-fl(i)).^2);
    fl(i)                   = f(b);
    fl_ii                   = b;

	[a b]                   = min((f-fr(i)).^2);
    fr(i)                   = f(b);
    fr_ii                   = b;
    A(i,fl_ii:(fr_ii-1))	= 1;
    fids(i,:)               = [fl_ii,fr_ii-1]; 
end

rnk         = sum(A, 2);
num_bands  	= find((rnk(2:end)>=rnk(1:(end-1))) & ...
                   (rnk(2:end)~=0)~=0, 1, 'last' )+1;
A           = A(1:num_bands, :);
cf          = cf(1:num_bands);
end


function x_stdft = stdft(x, N, K, N_fft)
% Compute the short-time hanning-windowed DFT of x. 
% Inputs:
%   N:          Frame size [samples].
%   K:          Overlap [samples].-, overlap K and DFT size
%   n_fft:      DFT size [samples]. 
% Outputs:
%   x_stdft:    Matrix with short-time DFT [-].
% The columns and rows of X_STDFT denote the frame-index and DFT-bin index, 
% respectively.

frames      = 1:K:(length(x)-N);
x_stdft     = zeros(length(frames), N_fft);

w           = hanning(N);
x           = x(:);

for i = 1:length(frames)
    ii              = frames(i):(frames(i)+N-1);
	x_stdft(i, :) 	= fft(x(ii).*w, N_fft);
end
end
end

