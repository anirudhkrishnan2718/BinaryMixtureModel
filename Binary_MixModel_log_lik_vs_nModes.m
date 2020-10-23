clc
clear all;

load('HMM_test_150_cells_binarized_spikes.mat', 'HMM_150_cells_binarized_spikes')
% sample_in = HMM_150_cells_binarized_spikes(1:100, :);
sample_in_unshuffled = HMM_150_cells_binarized_spikes(1:100, :);

nCells = size(sample_in_unshuffled, 1);
timeBins = size(sample_in_unshuffled, 2);



%% defining number of chunks and creating chunks. then assigning train and test sets

n_chunks = 3;
nTrials = 500;
timeBins_train = ((n_chunks - 1)/n_chunks) * timeBins;


% defining initial condition array
n_init_cond = 15;

% evenly splitting every Nth timebin into each chunk
sample_in_randomized = sample_in_unshuffled(:, randperm(timeBins));
chunk_1 = sample_in_randomized(:, 1:3:end);
chunk_2 = sample_in_randomized(:, 2:3:end);
chunk_3 = sample_in_randomized(:, 3:3:end);

% assigning train and test data sets
train_set_1 = [chunk_2 chunk_3];
test_set_1 = chunk_1;
train_set_2 = [chunk_3 chunk_1];
test_set_2 = chunk_2;
train_set_3 = [chunk_1 chunk_2];
test_set_3 = chunk_3;

%% running EM algo to find w_alpha and m_i_alpha for all N folds separately

n_modes_values = [25:1:30];
mode_vals = numel(n_modes_values);
avg_log_lik_overFolds = zeros(mode_vals, n_init_cond);


disp('function called')


for i = 1:1:mode_vals
    
    nModes_current = n_modes_values(1, i);
    w_IC_all = ones(1, nModes_current) ./ nModes_current;
    m_IC_array = 0.45 + (0.55 - 0.45) .* rand(nCells, nModes_current, n_init_cond);

    disp(['starting loop for nModes = ', num2str(nModes_current)])
    parfor j = 1:n_init_cond
        tic
        out_1 = nModes_current
        [out_2, out_3] = run_EM_algo_IC_init(train_set_1, nModes_current, nTrials, m_IC_array(:, :, j), w_IC_all);
        [out_4, out_5] = run_EM_algo_IC_init(train_set_2, nModes_current, nTrials, m_IC_array(:, :, j), w_IC_all);
        [out_6, out_7] = run_EM_algo_IC_init(train_set_3, nModes_current, nTrials, m_IC_array(:, :, j), w_IC_all);
        
        avg_log_lik_Fold1 = avg_log_lik_3Fold(out_2, out_3, test_set_1, nModes_current);
        avg_log_lik_Fold2 = avg_log_lik_3Fold(out_4, out_5, test_set_2, nModes_current);
        avg_log_lik_Fold3 = avg_log_lik_3Fold(out_6, out_7, test_set_3, nModes_current);
        
        avg_log_lik_overFolds(i, j) = (avg_log_lik_Fold1 + avg_log_lik_Fold2 + avg_log_lik_Fold3) / 3;
        toc
        datetime
    end
    
%     disp('parallelfor closing')
    
end

%% saving output per mode per distinct IC along with list of mode values used

disp('Saving results to mat fille')
save('CommonIC_15_Folds_3_EMIter_500_log_lik_vs_nModes_23_1_32.mat', 'n_modes_values', 'avg_log_lik_overFolds')
% loadedOutputs = load('CommonIC_9_Folds_3_EMIter_500_log_lik_vs_nModes.mat')



