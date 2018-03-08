%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Super Hexagon Statistical Analysis
% By: Matt Ford
% Date: February, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

% Load survival times
load time.mat

% Enable LaTeX formatting in figure text
set(0, 'defaulttextinterpreter', 'latex')

%% Look for drift in average performance
% Examine whether there are any systematic trends in game performance over
% time (this would violate the ergodic hypothesis). Find and plot the
% average performance in 10 evenly-spaced bins.

Nbins = 10;
binSize = length(time)/Nbins;

m = zeros(Nbins, 1);
mn = zeros(Nbins, 1);
mx = zeros(Nbins, 1);

for q = 1:Nbins
    
    i1 = floor((q-1)/Nbins * length(time)) + 1;
    i2 = ceil(q/Nbins      * length(time));
    
    time_q = time(i1:i2);
    time_q_sort = sort(time_q);
    m(q) = median(time_q);
    mn(q) = time_q_sort(4);
    mx(q) = time_q_sort(end-3);
    
end

% Plot results from each trial
figure(1)
plot((1:length(time))', time, 'k.',...
     'MarkerSize', 5, 'MarkerEdgeColor', [0.5 0.5 0.5]);
 
hold on
errorbar(linspace(0,length(time),Nbins)',m,m-mn,mx-m,'ko',...
         'MarkerFaceColor','k');

legend('all matches','median and 20% / 80% percentile')
     
axis([-5 155 0 max(time)+5])

xlabel('sample number', 'FontSize', 12);
ylabel('survival time [s]', 'FontSize', 12);

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition', [op(1:2) 450 350]);

% if ~exist('quantiles.png'); print(gcf,'quantiles.png','-dpng'); end

%% Look for correlation effects
% Is there a correlation between my performance on game i and game i+1? If
% so, this violate the hypothesis that the performance on subsequent games
% are independently distributed.

figure(2)
plot(time(1:end-1),time(2:end),'k.');
xlabel('survival time, match $n$','FontSize',12);
ylabel('survival time, match $n+1$','FontSize',12);

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 350]);

% if ~exist('next_test.png'); saveas(gcf,'next_test.png'); end

figure(3)
[acf, lags, bounds] = autocorr(time);

stem(lags, acf, 'filled', 'r', 'MarkerSize', 4); hold on
plot(lags, bounds(1)*ones(size(lags)), 'k:');
plot(lags, bounds(2)*ones(size(lags)), 'k:');
xlabel('lag', 'FontSize', 12);
ylabel('autocorrelation function', 'FontSize', 12);

op = get(gcf, 'OuterPosition');
set(gcf, 'OuterPosition', [op(1:2) 450 350]);

% if ~exist('autocorr.png'); saveas(gcf,'autocorr.png'); end

% Ljung-Box Q-Test
[h, p, Qstat, crit] = lbqtest(time, 'Lags', [1 2 3])

%% Plot histogram and estimate probability density function

[nh, th] = hist(time);
bar(th, nh, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k')

[nh_ks th_ks] = ksdensity(time, 'npoints', 100);
hold on
plot(th_ks,length(nh) * 120 * nh_ks * (th_ks(2) - th_ks(1)),'k');

xlabel('survival time','FontSize',12);
ylabel('counts','FontSize',12);

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 300]);

%% Initial model fitting
% Fit two statistical models to the data and make naive predictions.
% (1) The exponential model is based on the assumption that failure is a
%     Poisson process, and there is no time-dependence of my likelihood of
%     failure.
% (2) The Weibull model is an extension of the exponential distribution
%     which allows the likelihood of failure to change as a power law.

fexp = fitdist(time,'Exponential');
fwei = fitdist(time,'Weibull');
fwei2 = fitdist(time(1:120),'Weibull');

figure(4)
plot(sort(time),linspace(0,1,length(time)),'k.');

tt = linspace(0,max(time));
hold on;
plot(tt,1 - exp(-tt/fexp.mu),'k:','LineWidth',2)
plot(tt,1 - exp(-(tt/fwei.a).^fwei.b),'r--','LineWidth',2)

xlabel('survival time [s]','FontSize',12);
ylabel('cumulative probability','FontSize',12);
axis([0 max(time) 0 1.05])

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 350]);

legend('tests','exponential','Weibull','Location','best')

% if ~exist('cdf.png'); saveas(gcf,'cdf.png'); end

%% Failure probability with time

% Estimate the probability that I will survive to 60 seconds, based on both
% models
p60_wei = exp(-(60/fwei.a)^fwei.b);
p60_exp = exp(-(60/fexp.mu));

% Estimate the number of trials before my probability of success reaches
% 95%.
fconf = 0.95;
n_wei = log(1 - fconf)/log(1-p60_wei);
n_exp = log(1 - fconf)/log(1-p60_exp);

% Seconds to win with 95% probablity
T_wei = n_wei * (mean(time) + 1);
T_exp = n_exp * (mean(time) + 1);

fprintf('Number of years (Weibull): %.0f\n', (T_wei / 3600 / 24 / 365))
fprintf('Number of hours (Exp.)   : %3.2f\n', (T_exp / 3600))

%% Sensitivity analysis on weibull distribution fit
% Find predicted times to win based on upper and lower bounds of the
% confidence interval of the fitted Weibull model

a_max =   9.7258 + 0.429246;
b_max =  1.95849 + 0.120843;

a_min =   9.7258 - 0.429246;
b_min =  1.95849 - 0.120843;

p60_max = exp(-(60/a_max)^b_max);
p60_min = exp(-(60/a_min)^b_min);

n_max = log(1 - fconf)/(-p60_max);
n_min = log(1 - fconf)/(-p60_min);

T_max = n_max * (mean(time) + 1);
T_min = n_min * (mean(time) + 1);

(T_max / 3600 / 24 / 365) / 1e9
(T_min / 3600 / 24 / 365) / 1e9

figure(4)
plot(sort(time),linspace(0,1,length(time)),'k.');

tt = linspace(0,max(time));
hold on;
plot(tt,1 - exp(-(tt/fwei.a).^fwei.b),'r-','LineWidth',2)
plot(tt,1 - exp(-(tt/a_max).^b_max),'r:','LineWidth',1)
plot(tt,1 - exp(-(tt/a_min).^b_min),'r:','LineWidth',1)

xlabel('survival time [s]','FontSize',12);
ylabel('cumulative probability','FontSize',12);
axis([0 max(time) 0 1.05])

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 350]);

legend('tests','Weibull','confidence bounds','Location','best')

%% Sensitivity analysis
% Study variability of the predicted time to win based on which training
% data is used

[t s] = sort(rand(150,1));

w1 = fitdist(time(s(1:75)), 'Weibull');
w2 = fitdist(time(s(76:end)), 'Weibull');

p60_1 = exp(-(60/w1.a)^w1.b);
p60_2 = exp(-(60/w2.a)^w2.b);

n_1 = log(1 - fconf)/-p60_1;
n_2 = log(1 - fconf)/-p60_2;

T_1 = n_1 * (mean(time) + 1) / 3600 / 24 / 365
T_2 = n_2 * (mean(time) + 1) / 3600 / 24 / 365

%% Plot failure risk vs. time for Weibull and Exponential models

tt = linspace(0, max(time), 100);

plot(tt,(fwei.b/fwei.a)*(tt/fwei.a).^(fwei.b-1),'k','LineWidth',2);
hold on
plot(tt,ones(size(tt)) * 1/fexp.mu,'r--','LineWidth',2)

axis([0 max(time) 0 0.6])

xlabel('time [s]','FontSize',12);
ylabel('failure risk $f(t)$','FontSize',12);

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 300]);

legend('Weibull','exponential','Location','best')

%% Fit disjoint failure-risk model

T = sort(time);
F = linspace(0,0.99,length(time))';
T = T(4:end);
F = F(4:end);

% Fit an exponential model to the tail of the distribution, for t > 10
tc = 10;
ti = find(T > tc, 1);

pexp = polyfit(T(ti:end),log(1-F(ti:end)),1);
mu = -1/pexp(1);
Ic = pexp(2) - tc/mu;

plot(sort(time),linspace(0,1,length(time)),'k.');
hold on

% Plot exponential distribution fit
plot(tt(tt >= tc), 1-exp(Ic - (tt(tt >= tc) - tc)/mu),...
     'r', 'LineWidth', 2);
plot(tc, 1 - exp(Ic), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

xlabel('survival time [s]','FontSize',12);
ylabel('cumulative probability','FontSize',12);
axis([0 max(time) 0 1.05])

op = get(gcf,'OuterPosition');
set(gcf,'OuterPosition',[op(1:2) 450 350]);

legend('Tests', 'Disjoint model', 'Location', 'best')

p60 = exp(Ic - (60 - tc)*1/mu);

n_win = log(1 - fconf)/log(1-p60);
T_win = n_win * (mean(time) + 1);

fprintf('Days to ensure 95%% success: %.1f\n', T_win / 3600 / 24)
