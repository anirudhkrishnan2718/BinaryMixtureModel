function [out_unique_pop_response, out_emp_prob, out_model_prob, out_poisson_errorbar]= empirical_model_prob_unique_pop_response(w_a, m_i_a, test_set, nModes)
   
    %% find unique population responses and count number of occurences then find empirical prob
    tic
    timeBins_test = size(test_set, 2);
    
    uniqueFind_test_set = transpose(test_set);
    [C, ia, ic] = unique(uniqueFind_test_set, 'rows', 'sorted');
    freq_emp = (accumarray(ic, 1));
    poisson_error_bar = (freq_emp.^(0.5)) / timeBins_test;
    pop_emp_prob = freq_emp ./ timeBins_test;
    
    
    %% find model probability for each unique pop response
    
    test_set_unique = transpose(C);
    timeBins_test_unique = size(test_set_unique, 2);
    
    Q_test = zeros(nModes, timeBins_test_unique);
    P_mix_test = zeros(1, timeBins_test_unique);
    
    for u = 1:1:timeBins_test_unique
        for alpha = 1:1:nModes
            Q_test(alpha, u) = prod((m_i_a(:, alpha) .^ (test_set_unique(:, u))) .* ((1 - m_i_a(:, alpha)) .^ (1 - test_set_unique(:, u))));
        end
            P_mix_test(u) = sum(w_a.' .* Q_test(:, u));
    end
    
    
    %% return outputs
    
    out_unique_pop_response = test_set_unique;
    out_emp_prob = transpose(pop_emp_prob);
    out_model_prob = P_mix_test;
    out_poisson_errorbar = transpose(poisson_error_bar);
    
    disp('probability finding script done')
    toc
    
end