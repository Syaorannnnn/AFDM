function [X_hat] = Weighted_Hard_MRC_DFE(R, H_eff, N0, n_iter, M)
    
    % 功能: 真正的判决反馈均衡器 (引入了硬判决)
    % 输入:
    %   R       : 接收信号
    %   H_eff   : 等效信道矩阵
    %   N0      : 噪声功率
    %   n_iter  : 迭代次数
    %   M       : 调制阶数 (例如 4, 16, 64) 用于硬判决

    N_data = size(R,1);

    % 1. 预计算能量与稀疏索引 (保持不变)
    d = full(sum(abs(H_eff).^2, 1)).';
    col_rows = cell(N_data, 1);
    col_vals = cell(N_data, 1);
    for k = 1 : N_data
        [r, ~, v] = find(H_eff(:, k));
        col_rows{k} = r;
        col_vals{k} = v;
    end
    
    % 初始化
    x_soft = zeros(N_data, 1); % 软估计
    x_hard = zeros(N_data, 1); % 硬判决值
    Delta_y = R; 
    
    for n = 1 : n_iter
        for k = 1 : N_data
            rows = col_rows{k};
            if isempty(rows), continue; end
            h_vals = col_vals{k};
            
            % 加权 MRC (基于残差)
            % 注意：这里利用上一轮的残差来更新
            % 原始公式: g_k = h_k' * (y - sum_{j!=k} h_j x_j)
            % 变换为:   g_k = h_k' * (y - H*x_prev + h_k*x_prev_k) 
            %               = h_k' * Delta_y_prev + ||h_k||^2 * x_prev_k
            
            % 这里的 x_hard(k) 是上一轮(或本轮尚未更新)的硬判决值
            g_k = sum(conj(h_vals) .* Delta_y(rows)) + d(k) * x_hard(k);
            
            % LMMSE 软更新
            x_new_soft = g_k / (d(k) + N0);
            
            % --- 关键修改：硬判决 (Hard Decision) ---
            % 将软值映射到最近的星座点
            x_new_hard = qammod(qamdemod(x_new_soft, M, 'UnitAveragePower', true), ...
                                M, 'UnitAveragePower', true);
            
            % 计算增量 (用硬判决值的变化来更新残差)
            change = x_new_hard - x_hard(k);
            
            % 更新状态
            x_soft(k) = x_new_soft;
            x_hard(k) = x_new_hard;
            
            % 更新残差 (如果硬判决变了，才更新残差)
            if abs(change) > 1e-10
                Delta_y(rows) = Delta_y(rows) - h_vals * change;
            end
        end
    end
    
    % 输出：通常输出软值给解码器，或者直接输出硬值
    X_hat = x_soft; 
end