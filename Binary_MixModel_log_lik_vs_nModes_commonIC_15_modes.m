clc;

clear all;


%% load data from mat files

% input2 = load('CommonIC_15_Folds_3_EMIter_500_log_lik_vs_nModes_23_1_32.mat');
% input1 = load('CommonIC_15_Folds_3_EMIter_500_log_lik_vs_nModes_19_1_24.mat');

input1 = load('CommonIC_9_Folds_3_EMIter_500_log_lik_vs_nModes_25_5_50.mat');


%% using mean and standard deviation for plotting
% all_nModes = [input1.n_modes_values input2.n_modes_values];
% all_log_lik_raw_full = transpose([input1.avg_log_lik_overFolds ; input2.avg_log_lik_overFolds]);

all_nModes = [input1.n_modes_values];
all_log_lik_raw_full = transpose([input1.avg_log_lik_overFolds]);

num_nModes = numel(all_nModes);

%% using median and interpolated 68th percentile for plotting

percentile_threshold = 0.682689 % using mu + sigma as the errorbar threshold

avg_log_lik_medians = median(all_log_lik_raw_full, 1);
avg_log_lik_percentile_mu_sigma_upper = quantile(all_log_lik_raw_full, 0.5 + (0.5*percentile_threshold), 1);
avg_log_lik_percentile_mu_sigma_lower = quantile(all_log_lik_raw_full, 0.5 - (0.5*percentile_threshold), 1);

Y_ebar_pos = avg_log_lik_percentile_mu_sigma_upper - avg_log_lik_medians;
Y_ebar_neg = avg_log_lik_medians - avg_log_lik_percentile_mu_sigma_lower;
X_ebar_pos = zeros(1, num_nModes);
X_ebar_neg = zeros(1, num_nModes);

%% using raw data points for a scatter plot

all_nModes_scatter = repmat(all_nModes, 15, 1);
z1 = reshape(all_nModes_scatter, [], 1);
all_log_lik_raw_scatter = all_log_lik_raw_full;
z2 = reshape(all_log_lik_raw_scatter, [], 1);

%% plotting avg lok lik vs nModes using mean, stddev
% figure
% hold on
% grid on
% box on
% 
% avg_log_lik_chunks = mean(all_log_lik_raw_full, 1);
% avg_log_lik_stddev = std(all_log_lik_raw_full,0, 1);
% 
% 
% errorbar(all_nModes, avg_log_lik_chunks, avg_log_lik_stddev, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', 'black', 'MarkerEdgeColor', 'Red')
% % 
% % scatter(avg_log_lik(1, :), avg_log_lik(2, :), 'red', 'filled')
% % scatter(avg_log_lik(1, :), avg_log_lik(3, :), 'green', 'filled')
% % scatter(avg_log_lik(1, :), avg_log_lik(4, :), 'blue', 'filled')
% title('Average log likelihood using 3 fold cross validation and 9 distinct initial conditions as a function of number of modes', 'Interpreter','latex', 'FontSize', 20)
% xlabel('Number of modes', 'Interpreter','latex', 'FontSize', 14)
% ylabel('Average log likelihood from model', 'Interpreter','latex', 'FontSize', 14)
% % legend({'test chunk1', 'test chunk2', 'test chunk3'}, 'Interpreter','latex', 'FontSize', 12)
% legend({'mean and standard deviation for data points and errorbar'}) 


%% plotting avg lok lik vs nModes with median and mu +- sigma percentiles

f2 = figure;
hold on
grid on
box on

avg_log_lik_chunks = mean(all_log_lik_raw_full, 1);
avg_log_lik_stddev = std(all_log_lik_raw_full,0, 1);


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

set(f2,'PaperPositionMode','auto');         
set(f2,'PaperOrientation','landscape');
set(f2, 'PaperType', 'usletter');
print(f2, 'Ani_BinaryMixModel_logLik_vs_nModes', '-dpdf', '-fillpage');

%% plotting avg lok lik vs nModes raw data scatter plot
% figure
% hold on
% grid on
% box on
% 
% avg_log_lik_chunks = mean(all_log_lik_raw_full, 1);
% avg_log_lik_stddev = std(all_log_lik_raw_full,0, 1);
% 
% 
% scatter(z1, z2, 10, [0 0.5 0], 'filled')
% % 
% % scatter(avg_log_lik(1, :), avg_log_lik(2, :), 'red', 'filled')
% % scatter(avg_log_lik(1, :), avg_log_lik(3, :), 'green', 'filled')
% % scatter(avg_log_lik(1, :), avg_log_lik(4, :), 'blue', 'filled')
% title('Average log likelihood using 3 fold cross validation and 9 distinct initial conditions as a function of number of modes - raw data scatter', 'Interpreter','latex', 'FontSize', 20)
% xlabel('Number of modes', 'Interpreter','latex', 'FontSize', 14)
% ylabel('Average log likelihood from model', 'Interpreter','latex', 'FontSize', 14)
% % legend({'test chunk1', 'test chunk2', 'test chunk3'}, 'Interpreter','latex', 'FontSize', 12)
% legend({'raw data points for each distinct IC'})