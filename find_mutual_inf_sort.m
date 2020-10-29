function [out_mutual_inf_sorted, out_cell_index_identifier] = find_mutual_inf_sort(input_psth_data_set)

    all_cells_psth = input_psth_data_set;
    all_cells_psth_mean = mean(all_cells_psth, 2);
    all_cells_psth_normalized = all_cells_psth ./ all_cells_psth_mean;

    
    all_cells_m_i_terms = (log2(all_cells_psth_normalized) .* all_cells_psth);
%     ignore terms of this arraay which are NaN
    all_cells_m_i_terms(isnan(all_cells_m_i_terms)) = 0;
    
    all_cells_mutual_inf = sum(all_cells_m_i_terms, 2);
    [out_mutual_inf_sorted, out_cell_index_identifier] = sort(all_cells_mutual_inf, 'descend');
    
end
