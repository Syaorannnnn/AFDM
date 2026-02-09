function [H_eff_est, est_params_final] = Channel_Estimator(Y_rx, x_pilot, pilot_idx, N_data, CPP_len, c1, c2, P)
    % Channel_Estimator: AFDM 分数多普勒信道估计与矩阵重建
    %
    % 输入:
    %   Y_rx      : 接收到的 DAFT 域信号 (N_data x 1)
    %   x_pilot   : 发送的导频符号值 (标量)
    %   pilot_idx : 导频在 DAFT 帧中的位置索引 (1 ~ N_data)
    %   N_data    : 有效数据长度
    %   CPP_len   : CPP 长度 (用于确定时延搜索范围)
    %   c1, c2    : AFDM 参数
    %   P         : 需要估计的路径数量
    %
    % 输出:
    %   H_eff_est : 重建后的等效信道矩阵 (N_data x N_data)，可直接用于均衡
    %   est_params: 估计出的路径参数矩阵 [Delay, Doppler, ComplexGain]

    % 搜索配置
    range_l = 0 : CPP_len;              
    range_alpha = -2 : 2;
    steps = 1e-4;               
    grid_metrics = zeros(length(range_l), length(range_alpha));
    est_params_temp = zeros(P, 3); 
    
    % 粗搜索
    for i = 1 : length(range_l)
        for j = 1 : length(range_alpha)
            l = range_l(i);
            alpha = range_alpha(j);
            h_test = build_DAFT_response(N_data, c1, c2, l, alpha, pilot_idx - 1);
            grid_metrics(i, j) = abs(h_test' * Y_rx);
        end
    end
    
    temp_grid = grid_metrics;
    
    % 迭代提取路径
    for p = 1 : P
        [~, max_idx] = max(temp_grid(:));
        [r, c] = ind2sub(size(temp_grid), max_idx);
        
        hat_l = range_l(r);
        hat_alpha = range_alpha(c);
        
        % 精搜索
        best_metric = -inf;
        best_frac = 0;
        for frac = -0.5 : steps : 0.5
            nu_test = hat_alpha + frac;
            h_fine = build_DAFT_response(N_data, c1, c2, hat_l, nu_test, pilot_idx - 1);
            m = abs(h_fine' * Y_rx);
            if m > best_metric
                best_metric = m;
                best_frac = frac;
            end
        end
        hat_v = hat_alpha + best_frac;
        
        est_params_temp(p, :) = [hat_l, hat_v, 0];
        
        % 屏蔽邻域
        c_start = max(1, c - 2); 
        c_end   = min(size(temp_grid, 2), c + 2);
        temp_grid(r, c_start : c_end) = -inf; 
    end
    
    % LS 求解增益
    A_basis = zeros(N_data, P);
    for p = 1 : P
        l_p = est_params_temp(p, 1);
        f_p = est_params_temp(p, 2);
        A_basis(:, p) = build_DAFT_response(N_data, c1, c2, l_p, f_p, pilot_idx - 1);
    end
    
    A_matrix = A_basis * x_pilot;
    h_opt = pinv(A_matrix) * Y_rx;
    
    est_params_temp(:, 3) = h_opt;
    
    % 排序与相位去旋转
    % 按时延排序，方便对比
    [~, sort_idx] = sort(est_params_temp(:, 1)); 
    est_params = est_params_temp(sort_idx, :);
    
    % 去旋转 (Derotation)
    % 补偿 CPP 移除造成的相位偏移
    % N_data 归一化频率 -> 物理时间 CPP -> 相位
    derotation_phasor = exp(1j * 2 * pi * est_params(:, 2) * (CPP_len / N_data));
    
    est_params_final = est_params;
    est_params_final(:, 3) = est_params(:, 3) .* derotation_phasor;
    
    % 重建物理信道以及生成等效信道
    N_total = N_data + CPP_len;
    H_temp = Rebuild_LTVchannel(N_total, N_data, est_params_final);
    H_eff_est = Gen_Eff_Channel(H_temp, CPP_len, c1, c2);

end