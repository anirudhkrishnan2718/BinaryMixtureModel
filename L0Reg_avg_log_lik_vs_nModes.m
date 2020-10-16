clc;
clear all;

load('HMM_test_150_cells_binarized_spikes.mat', 'HMM_150_cells_binarized_spikes')
sample_in = HMM_150_cells_binarized_spikes(1:100, :);

nCells = size(sample_in, 1);
timeBins = size(sample_in, 2);



%% defining number of chunks and creating chunks. then assigning train and test sets

n_chunks = 3;

% evenly splitting every Nth timebin into each chunk
chunk_1 = sample_in(:, 1:3:end);
chunk_2 = sample_in(:, 2:3:end);
chunk_3 = sample_in(:, 3:3:end);

% assigning train and test data sets
train_set_1 = [chunk_2 chunk_3];
test_set_1 = chunk_1;
train_set_2 = [chunk_3 chunk_1];
test_set_2 = chunk_2;
train_set_3 = [chunk_1 chunk_2];
test_set_3 = chunk_3;



%% load results of EM algorithm as a struct
input_matFile_EM_results = load('threeFold_EM_Algo_500Iter_nModes_5_5_20.mat');

EMResults = input_matFile_EM_results.output_cellArray;
%% average log likelihood for each fold separately using L0 regularization term

n_modes_values = size(EMResults, 1);
nFolds = (size(EMResults, 2) - 1) / 2;
log_lik_array = zeros(n_modes_values, nFolds);
nModes_array = [];
etaVal = 0;

for i = 1:1:n_modes_values
    
    current_mode_val = EMResults{i, 1};
    nModes_array = [nModes_array; current_mode_val];
    current_w1 = EMResults{i, 2};
    current_m1 = EMResults{i, 3};
    current_w2 = EMResults{i, 4};
    current_m2 = EMResults{i, 5};
    current_w3 = EMResults{i, 6};
    current_m3 = EMResults{i, 7};
   
    avg_log_lik_L0_Fold1 = avg_log_lik_3Fold_L0Reg(current_w1, current_m1, test_set_1, etaVal, current_mode_val);
    avg_log_lik_L0_Fold2 = avg_log_lik_3Fold_L0Reg(current_w2, current_m2, test_set_2, etaVal, current_mode_val);
    avg_log_lik_L0_Fold3 = avg_log_lik_3Fold_L0Reg(current_w3, current_m3, test_set_3, etaVal, current_mode_val);
    
    log_lik_array(i, :) = [avg_log_lik_L0_Fold1 avg_log_lik_L0_Fold2 avg_log_lik_L0_Fold3];
    
end
    

%% plotting avg lok lik vs nModes with L0 reg term

figure
hold on
grid on
box on

avg_log_lik_chunks = mean(log_lik_array, 2);
avg_log_lik_stddev = std(log_lik_array,0, 2);


errorbar(nModes_array, avg_log_lik_chunks, avg_log_lik_stddev, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'Color', 'red')
% 
% scatter(avg_log_lik(1, :), avg_log_lik(2, :), 'red', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(3, :), 'green', 'filled')
% scatter(avg_log_lik(1, :), avg_log_lik(4, :), 'blue', 'filled')
title('Average log likelihood using 3 fold cross validation as a function of number of modes with L0 regularization $\eta = 0.0002$', 'Interpreter','latex', 'FontSize', 20)
xlabel('Number of modes', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average log likelihood from model', 'Interpreter','latex', 'FontSize', 14)
% legend({'test chunk1', 'test chunk2', 'test chunk3'}, 'Interpreter','latex', 'FontSize', 12)
legend({'test all chunks averaged'}) 