function h_vec = build_DAFT_response(N, c1, c2, l, v, tx_pilot_idx)
    % AFDM 理论响应向量 (用于估计阶段)
    p_indices = 0 : N-1; 
    q = tx_pilot_idx;    
    
    % 相位
    phase_term = (2*pi/N) * (N*c1*l^2 - q*l + N*c2*(q^2 - p_indices.^2));
    exp_phase = exp(1j * phase_term).';
    
    % 狄利克雷核
    theta = p_indices - q + v + 2*N*c1*l;
    numerator = exp(-1j * 2 * pi * theta) - 1;
    denominator = exp(-1j * 2 * pi/N * theta) - 1;
    
    F_term = zeros(N, 1);
    tol = 1e-6;
    idx_singularity = abs(denominator) < tol;
    F_term(~idx_singularity) = numerator(~idx_singularity) ./ denominator(~idx_singularity);
    F_term(idx_singularity) = N; 
    
    h_vec = (1/N) * exp_phase .* F_term;
    
end