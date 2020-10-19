clc;

clear all;


%% load data from mat files

input1 = load('CommonIC_9_Folds_3_EMIter_500_log_lik_vs_nModes.mat');
input2 = load("CommonIC_9_Folds_3_EMIter_500_log_lik_vs_nModes_25_5_50.mat");


%% using mean and standard deviation for plotting
all_nModes = [input1.n_modes_values input2.n_modes_values];
all_log_lik_raw = transpose([input1.avg_log_lik_overFolds; input2.avg_log_lik_overFolds]);

%% using median and interpolated 68th percentile for plotting

percentile_threshold = 0.682689 % using mu + sigma as the errorbar threshold

avg_log_lik_medians = median(all_log_lik_raw, 1);
avg_log_lik_percentile_mu_sigma_upper = quantile(all_log_lik_raw, 0.5 + (0.5*percentile_threshold), 1);
avg_log_lik_percentile_mu_sigma_lower = quantile(all_log_lik_raw, 0.5 - (0.5*percentile_threshold), 1);

Y_ebar_pos = avg_log_lik_percentile_mu_sigma_upper - avg_log_lik_medians;
Y_ebar_neg = avg_log_lik_medians - avg_log_lik_percentile_mu_sigma_lower;
X_ebar_pos = zeros(1, 10);
X_ebar_neg = zeros(1, 10);

%% plotting avg lok lik vs nModes using mean, stddev
figure
hold on
grid on
box on

avg_log_lik_chunks = mean(all_log_lik_raw, 1);
avg_log_lik_stddev = std(all_log_lik_raw,0, 1);


errorbar(all_nModes, avg_log_lik_chunks, avg_log_lik_stddev, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', 'black', 'MarkerEdgeColor', 'Red')
% 
% scatter(avg_log_lik(1, :), avg_log_lik(2, :), 'red', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(3, :), 'green', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(4, :), 'blue', 'filled')
title('Average log likelihood using 3 fold cross validation and 9 distinct initial conditions as a function of number of modes', 'Interpreter','latex', 'FontSize', 20)
xlabel('Number of modes', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average log likelihood from model', 'Interpreter','latex', 'FontSize', 14)
% legend({'test chunk1', 'test chunk2', 'test chunk3'}, 'Interpreter','latex', 'FontSize', 12)
legend({'mean and standard deviation for data points and errorbar'}) 


%% plotting avg lok lik vs nModes with median and mu +- sigma percentiles

figure
hold on
grid on
box on

avg_log_lik_chunks = mean(all_log_lik_raw, 1);
avg_log_lik_stddev = std(all_log_lik_raw,0, 1);


errorbar(all_nModes, avg_log_lik_medians, Y_ebar_neg, Y_ebar_pos, X_ebar_neg, X_ebar_pos, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', 'black', 'MarkerEdgeColor', 'blue')
% 
% scatter(avg_log_lik(1, :), avg_log_lik(2, :), 'red', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(3, :), 'green', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(4, :), 'blue', 'filled')
title('Average log likelihood using 3 fold cross validation and 9 distinct initial conditions as a function of number of modes', 'Interpreter','latex', 'FontSize', 20)
xlabel('Number of modes', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average log likelihood from model', 'Interpreter','latex', 'FontSize', 14)
% legend({'test chunk1', 'test chunk2', 'test chunk3'}, 'Interpreter','latex', 'FontSize', 12)
legend({'median and mu +- sigma percentiles for data points and errorbar'}) 