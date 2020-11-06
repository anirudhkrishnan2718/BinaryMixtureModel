clc
clear all;

load('HMM_test_300_cells_binarized_spikes.mat', 'HMM_300_cells_binarized_spikes')
sample_in_unshuffled = HMM_300_cells_binarized_spikes(1:200, :);

nCells = size(sample_in_unshuffled, 1);
timeBins = size(sample_in_unshuffled, 2);

% at frame 685 trial 1 begins
% at frame 70249 trial 50 ends

sample_in_unshuffled_firstHalf = sample_in_unshuffled(:, 685:70249);


%% defining number of chunks and creating chunks. then assigning train and test sets

nTrials = 500;
% number of EM algo iterations


% defining initial condition array
n_init_cond = 15;


train_set_1 = sample_in_unshuffled_firstHalf;
test_set_1 = sample_in_unshuffled_firstHalf;


%% running EM algo to find w_alpha and m_i_alpha for all N folds separately

n_modes_values = 35;
avg_log_lik_overFolds = zeros(1, n_init_cond);
% n_modes_values = [5:5:100]

disp('function called')
    
nModes_current = n_modes_values;
w_IC_all = ones(1, nModes_current) ./ nModes_current;
m_IC_array = 0.45 + (0.55 - 0.45) .* rand(nCells, nModes_current, n_init_cond);
w_dumpArray = {};
m_dumpArray = {};

datetime
disp(['starting loop for nModes = ', num2str(nModes_current)])

parfor j = 1:n_init_cond
    
    [out_2, out_3] = run_EM_algo_IC_init(train_set_1, nModes_current, nTrials, m_IC_array(:, :, j), w_IC_all);
    
    w_dumpArray{j, 1} = out_2;
    m_dumpArray{j, 1} = out_3;
    
    
    avg_log_lik_Fold1 = avg_log_lik_3Fold(out_2, out_3, test_set_1, nModes_current);
    avg_log_lik_overFolds(1, j) = avg_log_lik_Fold1;
    
    
end
datetime;
disp('parallelfor closing')
    

%% saving output per mode per distinct IC along with list of mode values used

save('CommonIC_15_Folds_1_EMIter_500_log_lik_vs_nModes_singleVal_35_firstHalf.mat', 'n_modes_values', 'avg_log_lik_overFolds', 'w_dumpArray', 'm_dumpArray')
save('TrainTestSets_RandomPermute_singleVal_35_firstHalf.mat', 'train_set_1', 'test_set_1')



