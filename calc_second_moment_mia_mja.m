function [out_emp_second_moment, out_model_second_moment, out_poi_ebar_second_moment] = calc_second_moment_mia_mja(w_a, m_i_a, test_set)
   
    model_second_moment = [];
    emp_second_moment = [];
    poi_ebar_second_moment = [];
    nCells = size(test_set, 1);
    timeBins_test = size(test_set, 2);
    
    for j = 1:1:nCells - 1
        for k = j+1:1:nCells
            
            newVal1 = sum(w_a .* m_i_a(k, :) .* m_i_a(j, :));
            newVal2 = sum(test_set(j, :) .* test_set(k, :)) ./ timeBins_test;
            newVal3 = sqrt(sum(test_set(j, :) .* test_set(k, :))) ./ timeBins_test;
            
            
            if sum(test_set(j, :) .* test_set(k, :)) > 0
                model_second_moment = [model_second_moment newVal1];
                emp_second_moment = [emp_second_moment newVal2];
                poi_ebar_second_moment = [poi_ebar_second_moment newVal3];
            else
                continue
                
            end
            
        end
    end
    
    out_emp_second_moment = emp_second_moment;
    out_model_second_moment = model_second_moment;
    out_poi_ebar_second_moment = poi_ebar_second_moment;
    
end