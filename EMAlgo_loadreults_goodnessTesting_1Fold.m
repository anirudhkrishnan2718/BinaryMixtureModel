clc
clear all;

load('HMM_test_300_cells_binarized_spikes.mat', 'HMM_300_cells_binarized_spikes')

sample_in_unshuffled = HMM_300_cells_binarized_spikes(1:200, :);

nCells = size(sample_in_unshuffled, 1);
timeBins = size(sample_in_unshuffled, 2);



%% defining number of chunks and creating chunks. then assigning train and test sets

nTrials = 500;
% number of EM algo iterations

n_chunks = 1;
n_modes_values = 35;

% defining initial condition array
n_init_cond = 15;


train_set_1 = sample_in_unshuffled;
test_set_1 = sample_in_unshuffled;
%% loading EM algo results and choosing a particular Initial condition to use results from

EMResults = load('CommonIC_15_Folds_1_EMIter_500_log_lik_vs_nModes_singleVal_35.mat');
median_log_lik = median(EMResults.avg_log_lik_overFolds);
median_IC_index = find(EMResults.avg_log_lik_overFolds == median_log_lik)

% manual override of which IC to use for computation
% median_IC_index = 7;

w_all = EMResults.w_dumpArray{median_IC_index, 1};
m_all = EMResults.m_dumpArray{median_IC_index, 1};

% fold by fold sorting w and m using w_alpha terms descending

w_input_unsorted_1 = w_all(1, :);
m_input_unsorted_1 = m_all(:, :, 1);

% sort w and m in descending order of weights
w_m_consolidated_1 = sortrows(transpose([w_input_unsorted_1; m_input_unsorted_1]), 1, 'descend');
w1 = transpose(w_m_consolidated_1(:, 1));
m1 = transpose(w_m_consolidated_1(:, 2:end));

%% comparing empirical vs model probability for each unique population response : calculation

[u_pop1, ep1, mp1, poi_ebar1] = empirical_model_prob_unique_pop_response(w1, m1, test_set_1, n_modes_values);



% leaving out the one count population responses
all_valid_emp_prob = [ep1(2:end)];
all_valid_model_prob = [mp1(2:end)];
all_poi_ebar = [poi_ebar1(2:end)];

% trimming points which have emp freq equal to 1
emp_model_prob_all = [all_valid_emp_prob ; all_valid_model_prob; all_poi_ebar];

threshold_emp_prob = 1.5 * (n_chunks/timeBins);
trimmed_data_emp_model_prob = [];
for b = 1:1:size(emp_model_prob_all, 2)
    if emp_model_prob_all(1, b) > threshold_emp_prob
       trimmed_data_emp_model_prob = [trimmed_data_emp_model_prob , emp_model_prob_all(:, b)];
    end
end

chi_sq_statistic_array = ((trimmed_data_emp_model_prob(1, :) - trimmed_data_emp_model_prob(2, :)).^2) ./ ((trimmed_data_emp_model_prob(3, :)).^2);
chi_sq_result = mean(chi_sq_statistic_array);

% disp({'Number of unique pop responses', size(all_valid_emp_prob, 2)});
% disp({'Number of pop responses with >1 frequency', size(trimmed_data_emp_model_prob, 2)});
% disp({'Chi squared statistic', chi_sq_result});
% disp({'(P_empirical - P_model)^2 / (Poisson_error_bar)^2 for each data point'})
% disp({'Poisson error bar = sqrt(frequency in empirical data) / total timeBins in empirical data'})
% disp({'Chi squared statistic is mean over all data points of this expression'})


%% comparing empirical vs model probability for each unique population response : plotting

f1 = figure;
hold on
grid on
box on


plot(trimmed_data_emp_model_prob(1, :), trimmed_data_emp_model_prob(2, :), 'Color', 'red', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue')
title({['Comparison of model probabilities to empirical probabilites'] ; ['for each unique population response (200 cells, 35 modes)'];[ '$\chi^{2} = $', num2str(chi_sq_result)]}, 'Interpreter','latex', 'FontSize', 20)
xlabel('Empirical Probability', 'Interpreter','latex', 'FontSize', 14)
ylabel('Model Probability', 'Interpreter','latex', 'FontSize', 14)

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')


refline(gca, 1, 0)

% create area plot using continuous error envelope to replace error bars

logspace_start = log10(threshold_emp_prob)
x_vals_env = logspace(logspace_start, -1, 1000);
y_vals_upper_env = x_vals_env + (sqrt(x_vals_env) ./ sqrt(timeBins / n_chunks));
y_vals_lower_env = x_vals_env - (sqrt(x_vals_env) ./ sqrt(timeBins / n_chunks));

x_vals_combined = [x_vals_env, fliplr(x_vals_env)];
inBetween = [y_vals_lower_env, fliplr(y_vals_upper_env)];
h1 = fill(x_vals_combined, inBetween, [0 0 0]);
set(h1, 'FaceAlpha', .25, 'EdgeAlpha', 0)



set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1, 'PaperType', 'usletter');
print(f1, 'Ani_BinaryMixModel_200Cell_35_Mode_Goodness_model_prob', '-dpdf', '-fillpage');






%% comparing empirical vs model first moment: calculation

[moment1_emp1, moment1_model1, p_ebar_fm_1] = calc_first_moment_m_i_alpha(w1, m1, test_set_1);

all_valid_emp_moment1 = [moment1_emp1];
all_valid_model_moment1 = [moment1_model1];
all_valid_moment1_poi_ebar = [p_ebar_fm_1];

chi_sq_statistic_array_fm = ((all_valid_emp_moment1 - all_valid_model_moment1).^2) ./ ((all_valid_moment1_poi_ebar).^2);
chi_sq_result_fm = mean(chi_sq_statistic_array_fm);


%% comparing empirical vs model first moment: plotting


f2 = figure;
hold on
grid on
box on

plot(all_valid_emp_moment1, all_valid_model_moment1, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'Color', 'red')
title({'Comparison of model to empirical first moments (200 cells, 35 modes) $ \sigma_{i} $', ['$ \chi^{2} =  $ ', num2str(chi_sq_result_fm)]}, 'Interpreter','latex', 'FontSize', 20)
xlabel('Empirical $m_{i \alpha}$', 'Interpreter','latex', 'FontSize', 14)
ylabel('Model $m_{i \alpha}$', 'Interpreter','latex', 'FontSize', 14)

refline(gca, 1, 0)

% create area plot using continuous error envelope to replace error bars

x_vals_env_fm = linspace(0, 0.035, 1000);
y_vals_upper_env_fm = x_vals_env_fm + (sqrt(x_vals_env_fm) ./ sqrt(timeBins / n_chunks));
y_vals_lower_env_fm = x_vals_env_fm - (sqrt(x_vals_env_fm) ./ sqrt(timeBins / n_chunks));

x_vals_combined_fm = [x_vals_env_fm, fliplr(x_vals_env_fm)];
inBetween_fm = [y_vals_lower_env_fm, fliplr(y_vals_upper_env_fm)];
h1_fm = fill(x_vals_combined_fm, inBetween_fm, [0 0 0]);
set(h1_fm, 'FaceAlpha', 0.1, 'EdgeAlpha', 0)

set(f2,'PaperPositionMode','auto');         
set(f2,'PaperOrientation','landscape');
set(f2, 'PaperType', 'usletter');
print(f2, 'Ani_BinaryMixModel_200Cell_35_Mode_Goodness_firstMoment', '-dpdf', '-fillpage');


%% comparing empirical vs model s moment: calculation

[moment2_emp1, moment2_model1, p_ebar_sm_1] = calc_second_moment_mia_mja(w1, m1, test_set_1);



all_valid_emp_moment2 = [moment2_emp1];
all_valid_model_moment2 = [moment2_model1];
all_valid_moment2_poi_ebar = [p_ebar_sm_1];

chi_sq_statistic_array_sm = ((all_valid_emp_moment2 - all_valid_model_moment2).^2) ./ ((all_valid_moment2_poi_ebar).^2);
chi_sq_result_sm = mean(chi_sq_statistic_array_sm);


%% comparing empirical vs model first moment: plotting


f3 = figure;
hold on
grid on
box on

plot(all_valid_emp_moment2, all_valid_model_moment2, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'Color', 'magenta');
title({'Comparison of model  to empirical second moments (200 cells, 35 modes) $ \sigma_{i}\sigma_{j}$', ['$ \chi^{2} =  $ ', num2str(chi_sq_result_sm)]}, 'Interpreter','latex', 'FontSize', 20)
xlabel('Empirical $ \sigma_{i}\sigma_{j}$', 'Interpreter','latex', 'FontSize', 14)
ylabel('Model $ \sigma_{i}\sigma_{j}$', 'Interpreter','latex', 'FontSize', 14)


% manual axis limits
xlim([0 1.8e-3])

refline(gca, 1, 0)


% create area plot using continuous error envelope to replace error bars

x_vals_env_sm = linspace(0, 6e-3, 1000);
y_vals_upper_env_sm = x_vals_env_sm + (sqrt(x_vals_env_sm) ./ sqrt(timeBins / n_chunks));
y_vals_lower_env_sm = x_vals_env_sm - (sqrt(x_vals_env_sm) ./ sqrt(timeBins / n_chunks));

x_vals_combined_sm = [x_vals_env_sm, fliplr(x_vals_env_sm)];
inBetween_sm = [y_vals_lower_env_sm, fliplr(y_vals_upper_env_sm)];
h1_sm = fill(x_vals_combined_sm, inBetween_sm, [0 0 0]);
set(h1_sm, 'FaceAlpha', .25, 'EdgeAlpha', 0)


set(f3,'PaperPositionMode','auto');         
set(f3,'PaperOrientation','landscape');
set(f3, 'PaperType', 'usletter');
print(f3, 'Ani_BinaryMixModel_200Cell_35_Mode_Goodness_secondMoment', '-dpdf', '-fillpage');


%% spikes per timebin histogram using synthetic data set : calculation

nTimeBins_test = (127500 / n_chunks);
syn_binary_spiketrain_test_1 = gen_synthetic_spiektrain_binary(w1, m1, nTimeBins_test);

spikes_per_timeBin_synth_1 = sum(syn_binary_spiketrain_test_1, 1);
[N1, x1] = histcounts(spikes_per_timeBin_synth_1);

spikes_per_timeBin_emp_1 = sum(test_set_1, 1);
[N2, x2] = histcounts(spikes_per_timeBin_emp_1);


%% spikes per timebin histogram using synthetic data set : plotting


f4 = figure;
hold on
grid on
box on

plot([1:1:numel(N1)] - 1, N1 ./ nTimeBins_test, 'r', 'LineStyle', ':', 'Marker', '.', 'MarkerSize', 14)
plot([1:1:numel(N2)] - 1, N2 ./ nTimeBins_test, 'b', 'LineStyle', '--', 'Marker', '.', 'MarkerSize', 14)
title('Population spike count distribution synthetic data vs test data set (200 cells, 35 modes)', 'Interpreter','latex', 'FontSize', 20)
xlabel('Spike count (k)', 'Interpreter','latex', 'FontSize', 14)
ylabel('Probability of k spikes $P_{k}$', 'Interpreter','latex', 'FontSize', 14)
legend({'synthetic data - trained using full data', 'empirical test - full data'}, 'Interpreter','latex', 'FontSize', 12)

set(gca, 'yscale', 'log')

set(f4,'PaperPositionMode','auto');         
set(f4,'PaperOrientation','landscape');
set(f4, 'PaperType', 'usletter');
print(f4, 'Ani_BinaryMixModel_200Cell_35_Mode_Goodness_SpikesPerTimebin', '-dpdf', '-fillpage');


