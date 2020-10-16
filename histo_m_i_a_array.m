clc

input_matFile_EM_results = load('threeFold_EM_Algo_500Iter_nModes_5_5_20.mat');
EMResults = input_matFile_EM_results.output_cellArray;

current_w1 = EMResults{4, 2};
current_m1 = EMResults{4, 3};
current_w2 = EMResults{4, 4};
current_m2 = EMResults{4, 5};
current_w3 = EMResults{4, 6};
current_m3 = EMResults{4, 7};


figure
hold on
grid on
box on

title('Histogram plot of m_i_a terms')
z1 = reshape(current_m1, 1, [])
histogram(z1)