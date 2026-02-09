function Channel_Estimation_Demo()
    clc; clear; close all;
    
    %% 1. 参数构造
    N = 256; 
    c1 = (2*(3+2)+1)/(2*N); % 简单设置
    c2 = 1/(N^2*2*pi);
    
    % 真实信道
    true_l = [0, 1, 2];
    true_nu = [2.45, -1.2, 0.3];
    true_h = [1.4, 1.1, 0.6] .* exp(1j*rand(1,3)); % 随机相位
    P = 3;
    
    % 构造接收信号 y (只含导频)
    x_pilot = 10;
    y = zeros(N, 1);
    for i = 1:P
        h_vec = build_daft_response_debug(N, c1, c2, true_l(i), true_nu(i), 0);
        y = y + true_h(i) * h_vec * x_pilot;
    end
    % 加点噪声
    y = y + 0.1 * (randn(N,1)+1j*randn(N,1));
    
    %% 2. 估计
    range_l = 0:3;
    range_alpha = -4:4;
    
    grid_metrics = zeros(length(range_l), length(range_alpha));
    
    % 粗搜索
    for i = 1:length(range_l)
        for j = 1:length(range_alpha)
            test_l = range_l(i);
            test_a = range_alpha(j);
            h_test = build_daft_response_debug(N, c1, c2, test_l, test_a, 0);
            grid_metrics(i,j) = abs(h_test' * y); % 必须取 abs!
        end
    end
    
    % 提取峰值
    est_params = zeros(P, 3);
    temp_grid = grid_metrics;
    
    fprintf('--- 开始提取 %d 条路径 ---\n', P);
    for p = 1:P
        [max_val, max_idx] = max(temp_grid(:));
        [r, c] = ind2sub(size(temp_grid), max_idx);
        
        hat_l = range_l(r);
        hat_a = range_alpha(c);
        
        fprintf('  路径 %d 粗搜索结果: Delay=%d, IntDoppler=%d (Corr=%.2f)\n', ...
                p, hat_l, hat_a, max_val);
        
        % 精搜索 (简化版)
        best_metric = -inf;
        best_frac = 0;
        for frac = -0.5:0.05:0.5
            h_fine = build_daft_response_debug(N, c1, c2, hat_l, hat_a+frac, 0);
            m = abs(h_fine' * y);
            if m > best_metric
                best_metric = m;
                best_frac = frac;
            end
        end
        hat_nu = hat_a + best_frac;
        
        % 估计 h
        h_final = build_daft_response_debug(N, c1, c2, hat_l, hat_nu, 0);
        hat_h = (h_final' * y) / (norm(h_final)^2 * x_pilot);
        
        est_params(p,:) = [hat_l, hat_nu, abs(hat_h)]; % 这里存模值方便对比
        
        % 确定屏蔽范围 (注意边界检查)
        idx_alpha_start = max(1, c - 1);
        idx_alpha_end   = min(size(temp_grid, 2), c + 1);
        
        % 将该时延下的邻近多普勒都设为极小值
        temp_grid(r, idx_alpha_start : idx_alpha_end) = -1;
    end
    
    %% 3. 结果对比
    fprintf('\n--- 最终结果 ---\n');
    disp('真实参数 (Delay, Doppler, Abs(Gain)):');
    disp([true_l; true_nu; abs(true_h)].');
    disp('估计参数 (Delay, Doppler, Abs(Gain)):');
    disp(est_params);
end

function h_vec = build_daft_response_debug(N, c1, c2, l, nu, q)
    p = (0:N-1).';
    phase = (2*pi/N)*(N*c1*l^2 - q*l + N*c2*(q^2 - p.^2));
    theta = p - q + nu + 2*N*c1*l;
    
    num = exp(-1j*2*pi*theta) - 1;
    den = exp(-1j*2*pi*theta/N) - 1;
    
    F = zeros(N,1);
    mask = abs(den)<1e-6;
    F(~mask) = num(~mask)./den(~mask);
    F(mask) = N;
    
    h_vec = (1/N) * exp(1j*phase) .* F;
end