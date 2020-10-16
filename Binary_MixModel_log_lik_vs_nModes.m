clc
clear all;

load('HMM_test_150_cells_binarized_spikes.mat', 'HMM_150_cells_binarized_spikes')
% sample_in = HMM_150_cells_binarized_spikes(1:100, :);
sample_in = HMM_150_cells_binarized_spikes(1:100, :);

nCells = size(sample_in, 1);
timeBins = size(sample_in, 2);



%% defining number of chunks and creating chunks. then assigning train and test sets

n_chunks = 3;
nTrials = 25;
timeBins_train = ((n_chunks - 1)/n_chunks) * timeBins;
output_cellArray = {};

% defining chunks

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

%% running EM algo to find w_alpha and m_i_alpha for all N folds separately

n_modes_values = [5]
% n_modes_values = [5:5:100]

disp('function called')

for i = 1:1:numel(n_modes_values)
    
    nModes = n_modes_values(1, i)
    output_cellArray{i, 1} = nModes;
    [output_cellArray{i, 2}, output_cellArray{i, 3}] = run_EM_algo_iffed(train_set_1, nModes, nTrials);
%     [output_cellArray{i, 4}, output_cellArray{i, 5}] = run_EM_algo_iffed(train_set_2, nModes, nTrials);
%     [output_cellArray{i, 6}, output_cellArray{i, 7}] = run_EM_algo_iffed(train_set_3, nModes, nTrials);
    
end

%% calculating log likelihood using EM algo output and held back test set incorporating cost term (L0 reg)


