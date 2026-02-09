function H = Rebuild_LTVchannel(N_total,N_data,params)
    % 重建物理信道 (线性移位模型)
    % params: [delay, doppler, gain]
    H = sparse(N_total, N_total);
    P = size(params, 1);
    
    doppler_cor = N_total / N_data;     % 修正多普勒参数,重建物理信道时应该使用的是相对于 N_total 的归一化多普勒参数

    for i = 1 : P
        l_i = round(params(i, 1));      % 时延必须是整数
        v_i = params(i, 2);             % 归一化多普勒 (相对于 N_data)
        h_i = params(i, 3);             % 复增益
        
        v_i = v_i * doppler_cor;        % 将多普勒参数从 N_data 归一化转换到 N_total 归一化
        rows = (1 + l_i) : N_total;
        cols = 1 : (N_total - l_i);
        Pi = sparse(rows, cols, ones(1, length(rows)), N_total, N_total);
        
        n_idx = (0 : N_total - 1).';
        
        D_diag = exp(-1j * 2 * pi * v_i * n_idx / N_total); 
        D_matrix = spdiags(D_diag, 0, N_total, N_total);
        
        H = H + h_i * D_matrix * Pi;
    end
end