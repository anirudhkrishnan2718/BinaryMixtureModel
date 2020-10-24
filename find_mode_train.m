function [out_dominant_mode] = find_mode_train(w_in, m_in, test_set_in)

    test_in = test_set_in;
    m_i_a = m_in;
    w_a = w_in;
    
    nCells = size(m_i_a, 1);
    nModes = size(m_i_a, 2);
    
    
    timeBins_test = size(test_in, 2);
    
    Q_test = zeros(nModes, timeBins_test);
    P_mix_indiv = zeros(nModes, timeBins_test);
    mode_train = zeros(1, timeBins_test);
    
    for u = 1:1:timeBins_test
        for alpha_test = 1:1:nModes
            Q_test(alpha_test, u) = prod((m_i_a(:, alpha_test) .^ (test_in(:, u))) .* ((1 - m_i_a(:, alpha_test)) .^ (1 - test_in(:, u))));
        end
        
        P_mix_indiv(:, u) = (w_a.') .* Q_test(:, u);
        [maxVal, maxIndex] = max(P_mix_indiv);
        
        out_dominant_mode = maxIndex;
        
    end
        
    
    
    
end