function [out_synthetic_binary_spiketrain] = gen_synthetic_spiektrain_binary(w_alpha, m_i_alpha, nTimeBins)
    tic
    nModes = size(m_i_alpha, 2);
    nCells = size(m_i_alpha, 1);
    syn_binary_spiketrain = -1 .* ones(nCells, nTimeBins);
    for t = 1:1:nTimeBins
        
        current_timeBin_modes = transpose(randsample([1:1:nModes], nCells, true, w_alpha));
        current_thresholds = -1*ones(nCells, 1);
        
        for n = 1:1:nCells
            current_thresholds(n, 1) = m_i_alpha(n, current_timeBin_modes(n, 1));
        end
        
        
        syn_binary_spiketrain(:, t) = (rand(nCells, 1) < current_thresholds);
        
    end
    
    out_synthetic_binary_spiketrain = syn_binary_spiketrain;
    toc
    disp('Synthetic spiketrains generated')
end