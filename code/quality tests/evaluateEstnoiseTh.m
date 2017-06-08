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
%%%%%%%%%%%%%%%%%%%%%%%%% th = 1.0 is the best %%%%%%%%%%%%%%%%%%%%%%%%%
targetSNR = [-20 -15 -10 -5 0 5 10 15 20];
targetSNR = targetSNR';

load 'snrStats_estnoise';
load 'pesqStats_estnoise';
load 'stoiStats_estnoise';

%% SNR results

pesqRatio_0pt7_avg = mean(pesqRatio_0pt7, 2);
pesqRatio_0pt8_avg = mean(pesqRatio_0pt8, 2);
pesqRatio_0pt9_avg = mean(pesqRatio_0pt9, 2);
pesqRatio_1pt0_avg = mean(pesqRatio_1pt0, 2);
pesqRatio_1pt1_avg = mean(pesqRatio_1pt1, 2);
pesqRatio_1pt2_avg = mean(pesqRatio_1pt2, 2);

%% segSNR results

segSNRratio_0pt7_avg = mean(segSNRratio_0pt7, 2);
segSNRratio_0pt8_avg = mean(segSNRratio_0pt8, 2);
segSNRratio_0pt9_avg = mean(segSNRratio_0pt9, 2);
segSNRratio_1pt0_avg = mean(segSNRratio_1pt0, 2);
segSNRratio_1pt1_avg = mean(segSNRratio_1pt1, 2);
segSNRratio_1pt2_avg = mean(segSNRratio_1pt2, 2);

%% STOI results

stoiRatio_0pt7_avg = mean(stoiRatio_0pt7, 2);
stoiRatio_0pt8_avg = mean(stoiRatio_0pt8, 2);
stoiRatio_0pt9_avg = mean(stoiRatio_0pt9, 2);
stoiRatio_1pt0_avg = mean(stoiRatio_1pt0, 2);
stoiRatio_1pt1_avg = mean(stoiRatio_1pt1, 2);
stoiRatio_1pt2_avg = mean(stoiRatio_1pt2, 2);

%%
segSNR = [segSNRratio_0pt7_avg, segSNRratio_0pt8_avg, segSNRratio_0pt9_avg, segSNRratio_1pt0_avg, segSNRratio_1pt1_avg, segSNRratio_1pt2_avg];

figure;
h = plot(targetSNR,segSNR,'linewidth',1.2);

title('\fontsize{19}Average segSNR ratios');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}segSNR ratio'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^';'d';'x'})
set(gca,'XTickLabel',aaaa,'fontsize',11.5)
legendCell = {'th=0.7', 'th=0.8', 'th=0.9', 'th=1.0', 'th=1.1', 'th=1.2'};
legend(legendCell,'FontSize',10.5)

%%
pesq = [pesqRatio_0pt7_avg, pesqRatio_0pt8_avg, pesqRatio_0pt9_avg, pesqRatio_1pt0_avg, pesqRatio_1pt1_avg, pesqRatio_1pt2_avg];

figure;
h = plot(targetSNR,pesq,'linewidth',1.2);

title('\fontsize{19}Average PESQ ratios');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}PESQ ratio'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^';'d';'x'})
set(gca,'XTickLabel',aaaa,'fontsize',11.5)
legendCell = {'th=0.7', 'th=0.8', 'th=0.9', 'th=1.0', 'th=1.1', 'th=1.2'};
legend(legendCell,'FontSize',10.5)

%%
stoi = [stoiRatio_0pt7_avg, stoiRatio_0pt8_avg, stoiRatio_0pt9_avg, stoiRatio_1pt0_avg, stoiRatio_1pt1_avg, stoiRatio_1pt2_avg];

figure;
h = plot(targetSNR,stoi,'linewidth',1.2);

title('\fontsize{19}Average STOI ratios');
xlabel('\fontsize{14}Global SNR of input noisy speech (dB)');
ylabel({'\fontsize{14}STOI ratio'});
aaaa = get(gca,'XTickLabel');
set(h,{'Marker'},{'s';'o';'*';'^';'d';'x'})
set(gca,'XTickLabel',aaaa,'fontsize',11.5)
legendCell = {'th=0.7', 'th=0.8', 'th=0.9', 'th=1.0', 'th=1.1', 'th=1.2'};
legend(legendCell,'FontSize',10.5)