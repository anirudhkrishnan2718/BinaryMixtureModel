function [out_emp_first_moment, out_model_first_moment, out_poi_ebar_first_moment] = calc_first_moment_m_i_alpha(w_a, m_i_a, test_set)
   
    nCells = size(test_set, 1)
    
    for k = 1:1:nCells
        
        model_first_moment(k) = sum(w_a .* m_i_a(k, :));

    end
    
    timeBins_test = size(test_set, 2);
    emp_first_moment = sum(test_set, 2) ./ timeBins_test;
    out_poi_ebar_first_moment = transpose((sum(test_set, 2) .^ 0.5)./ timeBins_test);
    
    
    out_model_first_moment = model_first_moment;
    out_emp_first_moment = transpose(emp_first_moment);
    
end