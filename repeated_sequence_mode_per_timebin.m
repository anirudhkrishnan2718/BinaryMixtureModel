clc;
clear all;

%% start and end of trimmed galvo signals
% 5064:7306089
% correspond to 127523 frames
% 75-25 split between 57ms and 58 ms frame time


%% isolate starting and ending frame indices
frames_per_sec = 17.5;
nModes = 35;
num_cells = 150;
num_reps = 90;
absolute_frame_positions = load('absolute_frame_positions_127523.mat');
frame_time_match = transpose(absolute_frame_positions.peak_locations);

start_end_times = load('frame_milliseconds_start_end_absolute.mat');
sorted_times = sort(start_end_times.Frame_start_end_milliseconds);

all_times = zeros(2, num_reps);
all_times_frames = zeros(2, num_reps);

for i = 1:1:num_reps
    all_times(1, i) = sorted_times(2*i - 1);
    all_times(2, i) = sorted_times(2*i);
end

gaps_millis = (all_times(2, :) - all_times(1, :));
gaps_frames = gaps_millis ./ 57;

%% find frame corresponding to start times by searching inside frame timing array

for k = 1:1:90
    target_start_time = all_times(1, k);
    target_end_time = all_times(2, k);
    
    all_times_frames(1, k) = find(frame_time_match == target_start_time);
    all_times_frames(2, k) = find(frame_time_match == target_end_time);
    
end
    
all_times_frames(3, :) = all_times_frames(2, :) - all_times_frames(1, :);

%% load results of EM algo for nModes single value (here n = 35)

EMResults = load('CommonIC_15_Folds_1_EMIter_500_log_lik_vs_nModes_singleVal_35.mat');
median_log_lik = median(EMResults.avg_log_lik_overFolds);
median_IC_index = find(EMResults.avg_log_lik_overFolds == median_log_lik);

% manual override of which IC to use for computation
median_IC_index = 7;

w_all = EMResults.w_dumpArray{median_IC_index, 1};
m_all = EMResults.m_dumpArray{median_IC_index, 1};


w_input_unsorted = w_all;
m_input_unsorted = m_all;

% sort w and m in descending order of weights
w_m_consolidated = sortrows(transpose([w_input_unsorted; m_input_unsorted]), 1, 'descend');
w_input = transpose(w_m_consolidated(:, 1));
m_input = transpose(w_m_consolidated(:, 2:end));


%% load cell arrays and calculate psth for all cells

all_raw_F = load('HMM_test_300_cells_binarized_spikes.mat');

all_cells = all_raw_F.HMM_300_cells_binarized_spikes(1:200, :);

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

%% plotting rannk ordered w_alpha terms for all 3 folds

f1 = figure;
hold on
grid on
box on


plot(transpose(sort(w_all, 2, 'descend')), 'Marker', '*', 'LineStyle', ':', 'Color', 'Red')

title('Rank ordered weights of all modes (200 cells 35 modes)', 'Interpreter','latex', 'FontSize', 20)
xlabel('mode index', 'Interpreter','latex', 'FontSize', 14)
ylabel('weight', 'Interpreter','latex', 'FontSize', 14)
legend({'test set 1 - full dataset'}, 'Interpreter','latex', 'FontSize', 14)
% saveas(gcf, 'Ani_BinaryMixModel_RankOrdered_ModeWeights.jpg')
set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1, 'PaperType', 'usletter');
print(f1, 'Ani_BinaryMixModel_RankOrdered_ModeWeights_200Cells_35Modes', '-dpdf', '-fillpage');


%% heatmap of mode per timebin across all trials of repeated stimulus

f2 = figure;
hold off
grid off
box on

% all modes simultaneously
colormap jet
h3 = imagesc(mode_train)
colorbar

axis([0 348 0 90])
title('Dominant mode at each timebin in repeated sequence trials (as colormap) (200 cells 35 modes)', 'Interpreter','latex', 'FontSize', 20)
xlabel('time bin in repeat sequence', 'Interpreter','latex', 'FontSize', 14)
ylabel('Trial index at each timebin', 'Interpreter','latex', 'FontSize', 14)
set(f2,'PaperPositionMode','auto');         
set(f2,'PaperOrientation','landscape');
set(f2, 'PaperType', 'usletter');
print(f2, 'Ani_BinaryMixModel_AllModes_200Cells_35Modes', '-dpdf', '-fillpage');

%% each individual mode as a subplot

f3 = figure;
sgtitle('Binary colormap for individual modes for each timebin in repeated sequence trials (200 cells 35 modes)', 'Interpreter','latex', 'FontSize', 20)
hold on
grid off
box on
colormap gray

axis([0 348 0 90])


for n = 1:1:nModes
    subplot(6, 6, n);
    imagesc(mode_train == n);
    set(gca, 'xtick', []);
    set(gca, 'xticklabel', []);
    set(gca, 'ytick', []);
    set(gca, 'yticklabel', []);
    
end

set(f3,'PaperPositionMode','auto');         
set(f3,'PaperOrientation','landscape');
set(f3, 'PaperType', 'usletter');
print(f3, 'Ani_BinaryMixModel_IndividualMode_Overview_200Cells_35Modes', '-dpdf', '-fillpage');

%% each mode in a new figure window 

for x = 1:1:nModes
    
    figure;
    ax1 = subplot(4, 4, [1:12]);
    colormap gray;
    imagesc(mode_train == x);
    set(ax1, 'xtick', []);
    set(ax1, 'xticklabel', []);
    set(ax1, 'ytick', []);
    set(ax1, 'yticklabel', []);
    
    ax2 = subplot(4, 4, [13:16]);
    plot(mode_emp_prob(x, :) ./ 0.057)
    
    xt = get(gca, 'XTick');                                 
    set(gca, 'XTick', xt, 'XTickLabel', xt/frames_per_sec)

    title(['Mode Index = ', num2str(x), ' and $w_{\alpha} = $', num2str(w_input(x))], 'Interpreter','latex', 'FontSize', 20)
    xlabel('time bin in repeat sequence', 'Interpreter','latex', 'FontSize', 14)
    ylabel('Mode Firing Rate', 'Interpreter','latex', 'FontSize', 14)
    
    filenameString = ['Ani_BinaryMixModel_IndividualMode_200Cells_35Modes_', num2str(x)]
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf, 'PaperType', 'usletter');
    print(gcf, filenameString, '-dpdf', '-fillpage');
    
    linkaxes([ax1 ax2], 'x');
end





%% plot each mode's empirical probability

f4 = figure;
hold on
grid on
box on
colormap jet
h2 = plot(transpose(mode_emp_prob ./ 0.057));

xt = get(gca, 'XTick');                                 
set(gca, 'XTick', xt, 'XTickLabel', xt/frames_per_sec)
    
title('Firing Rate of each mode across repeated sequence trials (200 cells 35 modes)', 'Interpreter','latex', 'FontSize', 20)
xlabel('time bin in repeat sequence', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average Firing Rate across all trials', 'Interpreter','latex', 'FontSize', 14)

set(f4,'PaperPositionMode','auto');         
set(f4,'PaperOrientation','landscape');
set(f4, 'PaperType', 'usletter');
print(f4, 'Ani_BinaryMixModel_AllModes_Probability_200Cells_35Modes', '-dpdf', '-fillpage');

hold off
