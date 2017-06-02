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

load 'snrStats_1';
load 'pesqStats_1';

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

%%
segSNR = [snrMaskLPC_avg, snrMDKFmask_avg, snrMDKF_avg, snrMMSE_avg, snrNoisy_avg];

figure;
h = plot(targetSNR,segSNR(:,[2,3,4,5]),'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'MDKF-IBM', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,segSNR(:,[1,3,4,5]),'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,segSNR,'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^';'d'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,segSNR(:,[1,2,3,4]),'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,segSNR(:,[2,3,4]),'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'MDKF-IBM', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,segSNR(:,[1,3,4]),'linewidth',1.2);

title('\fontsize{21}Average segSNR values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}segSNR (dB)'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)

%%
pesq = [pesqLPCmask_avg, pesqMDKFmask_avg, pesqMDKF_avg, pesqMMSE_avg, pesqNoisy_avg];

figure;
h = plot(targetSNR,pesq(:,[1,3,4,5]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^';'d'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq(:,[2,3,4,5]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^';'d'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'MDKF-IBM', 'MDKF', 'MMSE', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq(:,[1,3,4]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq(:,[2,3,4]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'MDKF-IBM', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq,'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^';'d'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF', 'MMSE', 'Noisy'};
% legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF', 'Noisy'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq(:,[1,2,3]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^'})
set(gca,'XTickLabel',aaaa,'fontsize',12)
legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF'};
legend(legendCell,'FontSize',11)

figure;
h = plot(targetSNR,pesq(:,[1,2,3,4]),'linewidth',1.2);

title('\fontsize{21}Average PESQ values');
xlabel('\fontsize{15}Global SNR of noisy speech (dB)');
ylabel({'\fontsize{15}PESQ'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'*';'^';'d'})
set(gca,'XTickLabel',aaaa,'fontsize',12)    
legendCell = {'LPC-enhanced', 'MDKF-IBM', 'MDKF', 'MMSE'};
legend(legendCell,'FontSize',11)
