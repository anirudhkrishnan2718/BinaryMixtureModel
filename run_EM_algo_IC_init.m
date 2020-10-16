function [w_output, m_i_a_output] = run_EM_algo_IC_init(input_sample_in, input_modes, nTrials, m_in, w_in)
%     tic
    
    trials = nTrials;
    sample_in = input_sample_in;
    timeBins = size(sample_in, 2);
    nCells = size(sample_in, 1);
    nModes = input_modes;
    w_a = w_in;    
    m_i_a = m_in;
    m_i_a_new = zeros(nCells, nModes);
    silent_state = zeros(nCells, 1);
    
    
    Q_array = zeros(nModes, timeBins);
    P_array = zeros(nModes, timeBins);

    
    for r = 1:1:trials
        
        
        m_i_a_silent_array = prod(1 - m_i_a, 1);
        
        %% expectation step
        for t = 1:1:timeBins
            for alpha = 1:1:nModes
                
                if isequal(sample_in(:, t), silent_state)
                    
%                     Q_array(alpha, t) = prod(1 - m_i_a(:, alpha));
                    Q_array(alpha, t) = m_i_a_silent_array(1, alpha);
                    
                else
                    
                    Q_array(alpha, t) = prod((m_i_a(:, alpha) .^ (sample_in(:, t))) .* ((1 - m_i_a(:, alpha)) .^ (1 - sample_in(:, t))));
                
                end
            end
            for alpha = 1:1:nModes
                
                P_array(alpha, t) = (w_a(alpha) .* Q_array(alpha, t)) / sum(w_a.' .* Q_array(:, t));
                
            end            
        end
        
        %% Maximization step
        w_a_new = transpose(sum(P_array, 2) / (timeBins + 1));
        for alpha = 1:1:nModes
            for i = 1:1:nCells
                m_i_a_new(i, alpha) = sum(sample_in(i, :) .* P_array(alpha, :)) / sum(P_array(alpha, :));
            end
        end
        
        w_a = w_a_new;
        m_i_a = m_i_a_new;
%         disp(['trial number ', num2str(r), ' for nModes = ', num2str(nModes)])
%         disp(['trial number ', num2str(r)])
    end
    
    

    
    
    w_output = w_a;
    m_i_a_output = m_i_a;
%     disp('EM_algo script done')
%     toc
    
end