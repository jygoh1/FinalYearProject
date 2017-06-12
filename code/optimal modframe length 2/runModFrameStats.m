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
load('pesqMDKF');
load('stoiMDKF');

pesqMDKF_avg = nanmean(pesqMDKF,2);
stoiMDKF_avg = nanmean(stoiMDKF,2);

Tw_slow = [16 20 24 32 64 128 256]*1e-3;

len = [];
for i = 1:length(Tw_slow)
    len{i} = num2str(Tw_slow(i)*1e3);
end


figure;
plot(pesqMDKF_avg,'linewidth',1.2);

set(gca, 'XTick', 1:length(Tw_slow), 'XTickLabel', len)
title('\fontsize{19}Average PESQ values');
xlabel('\fontsize{14}Modulation Frame Length (ms)');
ylabel({'\fontsize{14}PESQ'});
aaaa = get(gca,'XTickLabel');
set(gca,'XTickLabel',aaaa,'fontsize',12)


figure;
plot(stoiMDKF_avg,'linewidth',1.2);

set(gca, 'XTick', 1:length(Tw_slow), 'XTickLabel', len)
title('\fontsize{19}Average STOI values');
xlabel('\fontsize{14}Modulation Frame Length (ms)');
ylabel({'\fontsize{14}STOI'});
aaaa = get(gca,'XTickLabel');
set(gca,'XTickLabel',aaaa,'fontsize',12)
