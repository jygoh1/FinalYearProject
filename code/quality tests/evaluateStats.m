close all;
clearvars;
clc
beep off;

addpath(genpath('FYP'));
addpath(genpath('voicebox'));
addpath(genpath('pesqSTOI'));

%%
targetSNR = [-20 -15 -10 -5 0 5 10 15 20];
targetSNR = targetSNR';

% load 'snrStats_1';
% load 'pesqStats_1';
% load 'stoiStats_1';

load 'snrStats_2';
load 'pesqStats_2';
load 'stoiStats_2';

%% SNR results
snrMDKFmask_avg = mean(snrMDKFmask,2);
snrMaskLPC_avg = mean(snrMaskLPC,2);
snrMDKF_avg = mean(snrMDKF,2);
snrMMSE_avg = mean(snrMMSE,2);
snrNoisy_avg = mean(snrNoisy,2);

%% PESQ results
pesqLPCmask_avg = mean(pesqLPCmask, 2);
pesqMDKFmask_avg = mean(pesqMDKFmask, 2);
pesqMDKF_avg = mean(pesqMDKF, 2);
pesqNoisy_avg = mean(pesqNoisy, 2);
pesqMMSE_avg = mean(pesqMMSE, 2);

%% STOI results
stoiMDKFmask_avg = mean(stoiMDKFmask,2);
stoiLPCmask_avg = mean(stoiLPCmask,2);
stoiMDKF_avg = mean(stoiMDKF,2);
stoiMMSE_avg = mean(stoiMMSE,2);
stoiNoisy_avg = mean(stoiNoisy,2);

%%
segSNR = [snrMDKFmask_avg, snrMaskLPC_avg, snrMDKF_avg, snrMMSE_avg, snrNoisy_avg];

figure;
h = plot(targetSNR,segSNR(:,[1,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average segSNR values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'BMMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)

figure;
h = plot(targetSNR,segSNR(:,[2,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average segSNR values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)

%%
pesq = [pesqMDKFmask_avg, pesqLPCmask_avg, pesqMDKF_avg, pesqMMSE_avg, pesqNoisy_avg];

figure;
h = plot(targetSNR,pesq(:,[2,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average PESQ values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)

figure;
h = plot(targetSNR,pesq(:,[1,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average PESQ values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'BMMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)

%%
stoi = [stoiMDKFmask_avg, stoiLPCmask_avg, stoiMDKF_avg, stoiMMSE_avg, stoiNoisy_avg];

figure;
h = plot(targetSNR,stoi(:,[1,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average STOI values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'BMMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)

figure;
h = plot(targetSNR,stoi(:,[2,3,4,5]),'linewidth',1.2);

title('\fontsize{19}Average STOI values');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LMDKF', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',10.5)