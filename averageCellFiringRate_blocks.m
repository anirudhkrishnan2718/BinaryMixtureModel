clc;
clear all;

%% finding the number of frames in each kind of block

% full block = 1391;
% repeated_block = 346;
% random_block = 1045;

all_cells_spikes = load('AllCellsSpikeTrains_Steps1234_2370.mat');

binarized_spikes = all_cells_spikes.s_all_single_pre > 0;

nTimebins =  size(binarized_spikes, 2);
nCells =  size(binarized_spikes, 1);
num_reps = 90;
frames_per_sec = 17.5;

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
    


for m = 1:1:num_reps
   
    current_full_block = binarized_spikes(:, all_times_frames(2, m) - 1391 : all_times_frames(2, m));
    averageFiringRate(1, m) = find_average_firing_rate(current_full_block);
    
    current_random_block = binarized_spikes(:, all_times_frames(1, m) - 1045 : all_times_frames(1, m));
    averageFiringRate(2, m) = find_average_firing_rate(current_random_block);
    
    current_repeat_block = binarized_spikes(:, all_times_frames(1, m) : all_times_frames(1, m) + 346);
    averageFiringRate(3, m) = find_average_firing_rate(current_repeat_block);
    
end


%% plotting linked axes 3 plots

f1 = figure;
hold on

sgtitle('Average Cell Firing Rate in each block', 'Interpreter','latex', 'FontSize', 20)
ax1 = subplot(3, 3, [1:3]);
plot(averageFiringRate(1, :) .* frames_per_sec, 'Marker', '*', 'LineStyle', ':', 'Color', 'Red')
grid on
box on

% title(['Mode Index = ', num2str(x), ' and $w_{\alpha} = $', num2str(w_input(x))], 'Interpreter','latex', 'FontSize', 20)
% xlabel('time bin in repeat sequence', 'Interpreter','latex', 'FontSize', 14)
% ylabel('Mode Firing Rate', 'Interpreter','latex', 'FontSize', 14)
legend('Full block (80 seconds)', 'Interpreter','latex', 'FontSize', 12)

ax2 = subplot(3, 3, [4:6]);
plot(averageFiringRate(2, :) .* frames_per_sec, 'Marker', '*', 'LineStyle', ':', 'Color', 'Blue')
grid on
box on
% title(['Mode Index = ', num2str(x), ' and $w_{\alpha} = $', num2str(w_input(x))], 'Interpreter','latex', 'FontSize', 20)
% xlabel('time bin in repeat sequence', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average Cell Firing Rate (spikes per second)', 'Interpreter','latex', 'FontSize', 14)
legend('Random stimulus block (60 seconds)', 'Interpreter','latex', 'FontSize', 12)


ax3 = subplot(3, 3, [7:9]);
plot(averageFiringRate(3, :) .* frames_per_sec, 'Marker', '*', 'LineStyle', ':', 'Color', 'Black')
grid on
box on
% title(['Mode Index = ', num2str(x), ' and $w_{\alpha} = $', num2str(w_input(x))], 'Interpreter','latex', 'FontSize', 20)
xlabel('Block index', 'Interpreter','latex', 'FontSize', 14)
legend('Repeated stimulus block (20 seconds)', 'Interpreter','latex', 'FontSize', 12)


filenameString = ['Ani_BinaryMixModel_AverageFiringRate_per_block'];
set(f1,'PaperPositionMode','auto');         
set(f1,'PaperOrientation','landscape');
set(f1, 'PaperType', 'usletter');
print(f1, filenameString, '-dpdf', '-fillpage');

linkaxes([ax1 ax2 ax3], 'x');


%% plotting all 3 plots overlaid

f2 = figure;
hold on
grid on
box on

title('Average Cell Firing Rate in each block', 'Interpreter','latex', 'FontSize', 20);
plot(transpose(averageFiringRate) .* frames_per_sec, 'Marker', '*', 'LineStyle', ':');

xlabel('Stimulus block index', 'Interpreter','latex', 'FontSize', 14)
ylabel('Average Cell Firing Rate (spikes per second)', 'Interpreter','latex', 'FontSize', 14)
legend({'Full block (80 seconds)', 'Random stimulus block (60 seconds)', 'Repeated stimulus block (20 seconds)'}, 'Interpreter','latex', 'FontSize', 14)

filenameString = ['Ani_BinaryMixModel_AverageFiringRate_per_block_overlaid'];
set(f2,'PaperPositionMode','auto');         
set(f2,'PaperOrientation','landscape');
set(f2, 'PaperType', 'usletter');
print(f2, filenameString, '-dpdf', '-fillpage');