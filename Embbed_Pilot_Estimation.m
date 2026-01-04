function [theta_est] =  Embbed_Pilot_Estimation(y_pilot,N,c1,c2,P,max_Delay,max_Doppler)
    % 初始化网格搜索范围
    l_grid = 0 : max_Delay;
    alpha_grid = -ceil(max_Doppler) : ceil(max_Doppler);
    
    %% 粗估计
    % 假设分数部分为 0，遍历整数网格寻找能量峰值

    Metric_Map = zeros(length(l_grid), length(alpha_grid));
    
    for i_l = 1:length(l_grid)
        for i_a = 1:length(alpha_grid)
            l_trial = l_grid(i_l);
            nu_trial = alpha_grid(i_a); % 整数多普勒
            
            % 生成基向量 h(l, alpha)
            basis_vec = Gen_Pilot_Basis(N, c1, c2, l_trial, nu_trial);
            
            % 计算投影能量 |h^H * y|^2
            proj = basis_vec' * y_pilot;
            Metric_Map(i_l, i_a) = abs(proj)^2;
        end
    end
    
    % 寻找 P 个最大的峰值
    [~, sorted_indices] = sort(Metric_Map(:), 'descend');
    top_indices = sorted_indices(1:P);
    
    [l_idx, a_idx] = ind2sub(size(Metric_Map), top_indices);
    
    est_l_coarse = l_grid(l_idx);
    est_alpha_coarse = alpha_grid(a_idx);
    
    %% 细估计
    % 固定整数部分，搜索分数部分a[-0.5, 0.5]
    
    frac_grid = -0.5 : 0.05 : 0.5; % 搜索精度可调
    est_nu_final = zeros(P, 1);
    
    for i = 1:P
        best_metric = -inf;
        best_a = 0;
        
        l_curr = est_l_coarse(i);
        alpha_curr = est_alpha_coarse(i);
        
        for a_val = frac_grid
            nu_trial = alpha_curr + a_val;
            
            % 生成带分数多普勒的基向量
            basis_vec = Gen_Pilot_Basis(N, c1, c2, l_curr, nu_trial);
            
            % 计算 Metric (考虑归一化，虽然单位能量导频下分母近似为1)
            metric = abs(basis_vec' * y_pilot)^2 / (norm(basis_vec)^2);
            
            if metric > best_metric
                best_metric = metric;
                best_a = a_val;
            end
        end
        est_nu_final(i) = alpha_curr + best_a;
    end
    
    % 整理估计出的时延和多普勒
    est_delays = est_l_coarse(:);
    est_dopplers = est_nu_final(:);
    
    %% 第三步：增益求解 (Gain Estimation)
    % 由于分数多普勒破坏了正交性，需要解线性方程组
    % 对应论文 Eq (76): sum(h_j * correlation) = projection
    % 这等价于最小二乘问题: y = Phi * h
    
    Phi = zeros(N, P);
    for i = 1:P
        % 使用精细估计后的参数构建矩阵 Phi
        Phi(:, i) = Gen_Pilot_Basis(N, c1, c2, est_delays(i), est_dopplers(i));
    end
    
    % 求解 h_est = Phi \ y_pilot (最小二乘解)
    est_gains = Phi \ y_pilot;
    
    %% 输出结果
    theta_est.delays = est_delays;
    theta_est.dopplers = est_dopplers;
    theta_est.gains = est_gains;
end

%% 辅助函数：生成导频基向量
function h_vec = Gen_Pilot_Basis(N, c1, c2, l, nu)
    % 生成对应于参数 (l, nu) 的导频响应向量
    % 这是一个单位冲激导频通过信道 (l, nu) 后的 DAFT 域响应
    
    % 分解多普勒
    alpha = round(nu);
    frac = nu - alpha;
    
    % DAFT 域中的中心位置 loc (Eq 32)
    % 导频位于 index 0 (MATLAB index 1)
    % 注意：如果 c1 是基于 2N 的，这里公式要匹配
    loc = nu + 2 * N * c1 * l;
    
    h_vec = zeros(N, 1);
    
    % 确定能量泄漏范围 (kv)，分数多普勒通常取 1 或 2
    kv = 2; 
    
    % 构建基向量 (模拟 Eq 33-38 的 Sinc 扩散和相位旋转)
    for q_off = -kv : kv
        % 目标位置索引
        target_loc = round(loc) + q_off;
        q_idx = mod(target_loc, N); % 0-based index
        
        % 计算 Sinc 权重 (Dirichlet Kernel 近似)
        % frac部分实际上是 (loc - round(loc))
        dist = loc - target_loc; 
        weight = sinc(dist); % MATLAB sinc is sin(pi*x)/(pi*x)
        
        % 相位项 (Eq 33, 设 p=0)
        % phi = 2pi/N * (Nc1 l^2 - ql + Nc2 q^2)
        phase = (2*pi/N) * (N*c1*(l^2) - q_idx*l + N*c2*(q_idx^2));
        
        % 叠加到向量中 (1-based index)
        h_vec(q_idx + 1) = h_vec(q_idx + 1) + weight * exp(1j * phase);
    end

end