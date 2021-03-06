function [out_log_lik] = avg_log_lik_3Fold(w_in, m_in, test_set_in, nModes_in)

    test_in = test_set_in;
    m_i_a = m_in;
    w_a = w_in;
    nModes = nModes_in;
    
    timeBins_test = size(test_in, 2);
    
    Q_test = zeros(nModes, timeBins_test);
    P_mix_test = zeros(1, timeBins_test);
    
    for u = 1:1:timeBins_test
        for alpha_test = 1:1:nModes
            Q_test(alpha_test, u) = prod((m_i_a(:, alpha_test) .^ (test_in(:, u))) .* ((1 - m_i_a(:, alpha_test)) .^ (1 - test_in(:, u))));
        end
            P_mix_test(u) = sum(w_a.' .* Q_test(:, u));
    end
        
    out_log_lik = mean(log10(P_mix_test));
    
end