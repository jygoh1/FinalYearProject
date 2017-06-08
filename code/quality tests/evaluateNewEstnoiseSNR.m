close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

warning('off','all')
warning

%%
targetSNR = [-20 -15 -10 -5 0 5 10 15 20];
targetSNR = -targetSNR';

load 'estnoise_SNR';
load 'estnoiseIBM_SNR';

%% noise SNR results

noisepow_IBM_avg = mean(noisepow_IBM, 2);
noisepow_orig_avg = mean(noisepow_orig, 2);

%%
noiseSNR = [noisepow_orig_avg, noisepow_IBM_avg];

figure;
h = plot(targetSNR,noiseSNR,'linewidth',1.2);

hold on;
plot(targetSNR,targetSNR,'-.','linewidth',1.2,'Color',[0.2 0.8 0.2]);

title('\fontsize{19}Average estimated global SNR values');
xlabel('\fontsize{14}Actual global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}Estimated SNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*'})
set(h,{'Color'},{'r';'b'})
set(gca,'XTickLabel',aaaa,'fontsize',11.5)
legendCell = {'MMSE noise estimation','IBM-modified noise estimation'};
legend(legendCell,'FontSize',10.5)
