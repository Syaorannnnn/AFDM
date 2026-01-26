function [X_hat] = Weighted_MRC_DFE(R, H_eff, N0, n_iter)
    
    % 功能: 基于加权最大比合并的判决反馈均衡器 (Weighted MRC-DFE) 信号检测
    % 输入:
    %   R           : 去除了CPP的接收信号 (Nr x 1)
    %   H_eff       : 分数多普勒下的稀疏信道矩阵 (Nr x Nr)
    %   N0          : 噪声功率

    N_data = size(R,1);
    % H_trunc = H_eff(:, data_range); % 提取数据对应的信道列

    % 预计算每列能量
    d = full(sum(abs(H_eff).^2, 1)).';
    
    % 稀疏索引加速
    col_rows = cell(N_data, 1);
    col_vals = cell(N_data, 1);
    for k = 1 : N_data
        [r, ~, v] = find(H_eff(:, k));
        col_rows{k} = r;
        col_vals{k} = v;
    end
    
    x_curr = zeros(N_data, 1);
    x_prev = zeros(N_data, 1);
    Delta_y = R; 
    
    for n = 1 : n_iter
        for k = 1 : N_data
            rows = col_rows{k};
            if isempty(rows)
                continue; 
            end
            h_vals = col_vals{k};
            
            % 加权 MRC
            g_k = sum(conj(h_vals) .* Delta_y(rows)) + d(k) * x_prev(k);
            
            % LMMSE 更新
            x_curr(k) = g_k / (d(k) + N0);
            
            % 判决反馈
            change = x_curr(k) - x_prev(k);
            if abs(change) > 1e-8
                Delta_y(rows) = Delta_y(rows) - h_vals * change;
            end
        end
        x_prev = x_curr;
    end
    X_hat = x_curr;

end