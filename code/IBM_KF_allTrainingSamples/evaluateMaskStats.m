clear
clc
close all;

%%
load 'maskstats';

%%
figure
for i = 1:length(u_absent)
    mu = u_absent(i); sigma = var_absent(i);
    x = mu-(4*sigma):0.01:mu+(4*sigma);
    y = normpdf(x,mu,sigma);
    plot(x,y,'b')
    hold on;
end

figure
for i = 1:length(u_absent)
    mu = u_present(i); sigma = var_present(i);
    x = mu-(4*sigma):0.01:mu+(4*sigma);
    y = normpdf(x,mu,sigma);
    plot(x,y,'r')
    hold on;
end

figure
i = 93;
mu = u_absent(i) 
sigma = var_absent(i);
x = mu-(4*sigma):0.01:mu+(4*sigma);
y = normpdf(x,mu,sigma);
plot(x,y,'b','Linewidth',1.1)
hold on;
mu = u_present(i) 
sigma = var_present(i);
x = mu-(4*sigma):0.01:mu+(4*sigma);
y = normpdf(x,mu,sigma);
plot(x,y,'r','Linewidth',1.1)
hold on;
legend('Noise-dominant PDF','Speech-dominant PDF');
title('Probability distributions for noise-dominant and speech-dominant regions');
