clc;
clear all;

%% defining constants

nModes = 25;
nCells = 100;
num_reps = 90;


%% load mat files 

cell_psth_vals_struct = load('all2370Cells_repeat_segment_psth.mat');
all_times_frames_struct = load('repeat_segment_start_end_frames.mat');

cell_psth_vals = cell_psth_vals_struct.all_cells_psth;
all_times_frames = all_times_frames_struct.all_times_frames;


%% load results of EM algo for nModes single value (here n = 25) includes sorting modes by w_alpha magnitude

EMResults = load('CommonIC_15_Folds_3_EMIter_500_log_lik_vs_nModes_singleVal_25.mat');
median_log_lik = median(EMResults.avg_log_lik_overFolds);
median_IC_index = find(EMResults.avg_log_lik_overFolds == median_log_lik);

% manual override of which IC to use for computation
% median_IC_index = 7;

w_all = EMResults.w_dumpArray{median_IC_index, 1};
m_all = EMResults.m_dumpArray{median_IC_index, 1};

fold_index = 3;
w_input_unsorted = w_all(fold_index, :);
m_input_unsorted = m_all(:, :, fold_index);

% sort w and m in descending order of weights
w_m_consolidated = sortrows(transpose([w_input_unsorted; m_input_unsorted]), 1, 'descend');
w_input = transpose(w_m_consolidated(:, 1));
m_input = transpose(w_m_consolidated(:, 2:end));

%% load cell arrays and calculate mode sequence for all repeated segments

all_raw_F = load('HMM_test_150_cells_binarized_spikes.mat');

all_cells = all_raw_F.HMM_150_cells_binarized_spikes(1:100, :);

for m = 1:1:num_reps

        current_start = all_times_frames(1, m);
        z_repeated_new_trial = all_cells(:, current_start : current_start + 346);
        mode_train(m, :) = find_mode_train(w_input, m_input, z_repeated_new_trial);

end


%% find empirical probability of each mode per timebin across trials

for i = 1:1:nModes
    
    currentModeCheck = mode_train == i;
    current_emp_prob = sum(currentModeCheck, 1) ./ num_reps;
    mode_emp_prob(i, :) = current_emp_prob;

end

%% scatter mutual_inf vs average activity for cells

[z1, z2] = find_mutual_inf_sort(cell_psth_vals);
cell_psth_means = mean(cell_psth_vals, 2);
z3 = mean(all_cells, 2);
for m = 1:1:numel(z2)
   cell_psth_mean_sorted_mutual_inf(m, 1) = cell_psth_means(z2(m));
end

[z4, z5] = find_mutual_inf_sort(mode_emp_prob);
mode_emp_prob_means = mean(mode_emp_prob, 2);
z6 = mean(mode_emp_prob, 2);
for n = 1:1:numel(z5)
   mode_emp_prob_means_sorted_mutual_inf(n, 1) = mode_emp_prob_means(z5(n));
end

%% plotting for cells

f1 = figure;
hold on
grid on
box on


scatter(cell_psth_mean_sorted_mutual_inf(1:100)./0.057, z1(1:100), 20, 'Red', 'filled')

title('Mutual information vs average activity for cells', 'Interpreter','latex', 'FontSize', 20)
xlabel('average activity (spikes per second)', 'Interpreter','latex', 'FontSize', 14)
ylabel('Mutual information (bits per second)', 'Interpreter','latex', 'FontSize', 14)
% legend({'test set 1', 'test set 2', 'test set 3'})
% saveas(gcf, 'Ani_BinaryMixModel_RankOrdered_ModeWeights.jpg')
set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1, 'PaperType', 'usletter');
print(f1, 'Ani_BinaryMixModel_Scatter_mutualInf_vs_AvgActivity_Cells', '-dpdf', '-fillpage');

%% plotting for modes


f1 = figure;
hold on
grid on
box on

scatter(mode_emp_prob_means_sorted_mutual_inf, z4, 20, 'Blue', 'filled')

title('Mutual information vs average activity for modes', 'Interpreter','latex', 'FontSize', 20)
xlabel('average empirical probability', 'Interpreter','latex', 'FontSize', 14)
ylabel('Mutual information (bits per seond)', 'Interpreter','latex', 'FontSize', 14)
% legend({'test set 1', 'test set 2', 'test set 3'})
% saveas(gcf, 'Ani_BinaryMixModel_RankOrdered_ModeWeights.jpg')
set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1, 'PaperType', 'usletter');
print(f1, 'Ani_BinaryMixModel_Scatter_mutualInf_vs_AvgActivity_Modes', '-dpdf', '-fillpage');