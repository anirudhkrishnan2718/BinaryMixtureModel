function out_silentStateFraction = find_Silent_State_Fraction(input_array)
   
    spikes_block = input_array;
    nCells = size(spikes_block, 1);
    nTimeBins = size(spikes_block, 2);
    
    spikes_array = sum(spikes_block, 1);
    silent_state_num = nnz(~spikes_array);
    
    % transposeing and using ismember with rows argument
    out_silentStateFraction = silent_state_num / numel(spikes_array);   
    
end