function Heff = diag_reshape(y, pilot_in, gate, N, max_Delay, max_Doppler_alpha, c_1, c_2, pilot_loc, k_v)
% diag_reshape_v2: 基于聚类的鲁棒重构算法

%% 0. 预处理与参数准备
if length(pilot_in) > 1
    pilot_val = pilot_in(1);
else
    pilot_val = pilot_in;
end

% 动态门限 (关键：防止固定门限在低SNR下失效或高SNR下漏检)
% 如果 gate 输入为 0，则自动计算；否则使用输入值
if gate == 0
    max_val = max(abs(y));
    real_gate = max_val * 0.1; % 设定为最大峰值的 10%
else
    real_gate = gate;
end

% 搜索区域 (含 mod 防越界)
guardWidth = (max_Delay + 1) * (2 * (max_Doppler_alpha + k_v) + 1) - 1;
raw_search_idx = (max_Doppler_alpha + k_v) : (guardWidth + max_Doppler_alpha + k_v);
Pilot_area = mod(raw_search_idx - 1, N) + 1; 

%% 1. 寻峰与聚类 (Clustering)
% 不再对所有过门限的点一视同仁，而是先找“物理路径中心”

Heff = zeros(N,N);
y_search = zeros(size(y));
y_search(Pilot_area) = abs(y(Pilot_area)); % 只在合法区域搜索

% 寻找局部最大值作为路径中心 (Path Centers)
[pks, locs] = findpeaks(y_search, 'MinPeakHeight', real_gate, 'MinPeakDistance', 2*k_v);

% 如果没有 findpeaks 工具箱，可以用简单的循环替代，但建议保留

%% 2. 对每个簇进行统一重构
for i = 1:length(locs)
    center_idx = locs(i); % 簇中心的索引 (1-based)
    
    % --- A. 估计该簇的统一时延 (Robust Delay Estimation) ---
    % 仅根据能量最强的中心点来判断时延，最准确
    m = mod(center_idx - 1 - (2*max_Delay*max_Doppler_alpha + 2*max_Doppler_alpha + max_Delay), N); % 映射回逻辑位置
    
    d_est = 0;
    min_dist = inf;
    
    for d = 0:max_Delay
        % 理论中心位置
        theory_loc = mod(-2*N*c_1*d, N);
        
        % 环形距离计算
        dist = abs(m - theory_loc);
        if dist > N/2, dist = N - dist; end
        
        if dist < min_dist
            min_dist = dist;
            d_est = d;
        end
    end
    
    % 安全检查：如果距离太远，可能是虚警，跳过
    if min_dist > (max_Doppler_alpha + k_v + 2)
        continue;
    end
    
    % --- B. 扩展并重构该簇的所有带状线 ---
    % 强制取左右 k_v 个点，全部认为是属于时延 d_est 的分数多普勒分量
    cluster_range = -k_v : k_v;
    
    for offset = cluster_range
        % 当前要处理的“种子”位置
        curr_seed_idx = mod(center_idx + offset - 1, N) + 1;
        
        % 获取种子值 (Gain * Phase_start)
        h_seed = y(curr_seed_idx) / pilot_val;
        
        % 填入种子到矩阵第一列 (pilot_loc)
        Heff(curr_seed_idx, pilot_loc) = h_seed;
        
        % --- C. 向量化相位生成 (替代循环递推，速度快且无累积误差) ---
        % 我们需要计算对角线上所有点的值。
        % AFDM 沿路径的相位变化规律是确定的。
        
        % 生成 1 到 N-1 的步长向量 (代表从 pilot_loc 往后走了多少步)
        k_steps = (1:N-1).'; 
        
        % 计算行索引向量 p_vec 和 列索引向量 q_vec
        % 种子位置 (0-based)
        p_seed = curr_seed_idx - 1;
        q_seed = pilot_loc - 1;
        
        % 轨迹上的所有点 (0-based)
        p_vec = mod(p_seed + k_steps, N);
        q_vec = mod(q_seed + k_steps, N);
        
        % 直接计算相位偏移 (Direct Phase Calculation)
        % 公式推导自递推关系：
        % Phase_diff(k) = (2pi/N) * [ -d * k + N*c2 * ( (q^2 - p^2) - (q_seed^2 - p_seed^2) ) ]
        % 注意：这里 q - q_seed = k (如果不考虑循环), 
        % 考虑到 mod N，我们直接用坐标计算最稳妥
        
        term1 = -d_est * (q_vec - q_seed); % 线性项 (可能需要处理循环跳变，但 c2 项主导)
        
        % 更稳健的相位重建方法：利用 chirps
        % H(p,q) ∝ exp( j*2pi * (c1*l^2 - q*l/N + c2*q^2 - c2*p^2) )
        % 因此 H(p,q) / H(p_seed, q_seed) = ...
        
        % 为避免复杂的解包裹，我们还是用递推公式的向量化版本
        % 每一对 (p,q) 和前一对的差值
        % phase_step = (q_prev - q)*d + N*c2*(q^2 + p_prev^2 - p^2 - q_prev^2)
        % 这是一个随 k 变化的数列，我们可以累加它
        
        % 使用循环填充这一条对角线 (为了稳健性，不用纯向量化)
        curr_val = h_seed;
        prev_p = p_seed;
        prev_q = q_seed;
        
        for step = 1:N-1
            % 当前坐标
            curr_p = mod(prev_p + 1, N);
            curr_q = mod(prev_q + 1, N);
            
            % 1-based 索引用于矩阵赋值
            mat_p = curr_p + 1;
            mat_q = curr_q + 1;
            
            % 计算相位旋转 (Eq. 30)
            % 注意：此处 d_est 是该簇统一且唯一的时延！这是改进的关键！
            phase_term = (prev_q - curr_q)*d_est + N*c_2*(curr_q^2 + prev_p^2 - curr_p^2 - prev_q^2);
            tao = exp(1i * 2 * pi * phase_term / N);
            
            % 更新值
            curr_val = curr_val * tao;
            Heff(mat_p, mat_q) = curr_val;
            
            % 迭代
            prev_p = curr_p;
            prev_q = curr_q;
        end
    end
end

end