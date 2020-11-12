function [out_avg_firing_rate] = find_average_firing_rate(in_test_set)
    
    spikes_block = in_test_set;
    
    nCells = size(spikes_block, 1);
    nTimeBins = size(spikes_block, 2);
    
    % mean2 finds mean of all elements in array. mathematically the same as
    % finding average of each row and then average of all those averages
    
    out_avg_firing_rate = mean2(spikes_block);
    
    
end
