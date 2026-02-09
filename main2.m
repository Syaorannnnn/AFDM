clc; clear all; close all;
%% 基础参数
    bps = 2;    % 比特率
    M = 2^bps;  % 星座图点数
    N_data = 2^8;  % 有效数据数
    SNRdB = 0 : 5 : 20;  % 信噪比
    SNR_idx_max = length(SNRdB);

    P = 3; % 多径数
    N0 = 10.^(-SNRdB ./ 10);

    target_error_bits = 1e3;
    max_sym = 1e6;

    maxNormDoppler = 2; % 4GHz 载波下的 540 km/h
    % 给出固定的归一化时延，初始化归一化多普勒、信道增益
    pathDelays = [0 1 2];
    pathDopplers = zeros(P);
    pathGains = zeros(P);

    maxDelay = max(pathDelays);

    CPP_lenth = maxDelay;

    N_total = N_data + CPP_lenth; % 所需子载波个数 Data + Overhead
    data_range = (1 + CPP_lenth) : N_total;

    % 设计c1、c2参数、设计帧结构 
    k_v = 2;    % 用于对抗分数多普勒，为不同径带来足够高的分离度
    c1 = (2 * (maxNormDoppler + k_v) + 1) / (2 * N_data);           % 满足正交条件的c1，使多径各路径不发生交叠
    c2 = 1 / (N_data^2 * 2 * pi);
    
    if (2 * (maxNormDoppler + k_v) * (maxDelay + 1)) + maxDelay > N_data      %必须满足正交条件
        fprintf("子载波不满足正交！\n");
    end

    Q_calc = (maxDelay + 1) * (2 * (maxNormDoppler + k_v) + 1) - 1;
    Q = ceil(Q_calc);

    % 导频设计
    pilot_idx = 1;
    pilot = 10 + 0j;

    if pilot_idx + 2 * Q >= N_data
        fprintf("N_data过小,放不下双边间隔! \n");
    end
%% 蒙特卡洛仿真
    %% 多帧模拟
        BER = zeros(1,SNR_idx_max);    % 初始化BER数组
        %% SNR步进
        for i = 1 : SNR_idx_max
            % 初始化
            total_errs = 0;
            total_bits = 0;
            sym_count = 0;
            fprintf("当前正在仿真 %d SNRdB ... \n", SNRdB(i));
            while(total_errs <= target_error_bits && sym_count < max_sym)
                % Jakes分布
                % P个随机相位theta
                theta = (rand(1, P) * 2 * pi) - pi;
                pathDopplers = maxNormDoppler * cos(theta);
                
                % 瑞利信道
                pathGains = (randn(1, P) + 1j * randn(1, P)) / sqrt(2 * P);

                max_Doppler = max(abs(pathDopplers)); 
                max_Doppler_alpha = round(max_Doppler);                 % 对多普勒取整
                max_Doppler_a = max_Doppler - max_Doppler_alpha;        % 多普勒的小数部分，范围[-0.5,0.5]
                    
                    % 有效数据个数
                    N_active = N_data - 2 * Q - 1;  % 减1是因为还要给导频留位置
                    active_idx = (Q + 2) : (Q + 1 + N_active); % 有效数据位置索引

                    if max(active_idx) > (N_data - Q)
                        error('索引分配错误');
                    end

                    H = LTVchannel(N_total,pathDelays,pathDopplers,pathGains);                % 构造LTV信道矩阵
                    % H_eff_ref = Gen_Eff_Channel(H,CPP_lenth,c1,c2);                         % AFDM的有效信道矩阵
                    % H_eff = Gen_Eff_Channel(H,CPP_lenth,0,0);                         % OFDM的有效信道矩阵

                    % 每生成一个符号，计数器+1
                    sym_count = sym_count + 1;
    
                %% 发送端
                    % 生成原始数据
                    X0_active = randi([0 M-1], N_active, 1);
                    % 映射成星座点
                    X_QAM_active = qammod(X0_active, M, 'UnitAveragePower', true);
                    % 组建完整的DAFT域符号向量 (加入ZP & pilot)
                    X_QAM = zeros(N_data,1);
                    X_QAM(pilot_idx) = pilot;
                    X_QAM(active_idx) = X_QAM_active;
                    % IDAFT(c1 = c2 = 0时IDAFT退化成为IDFT)
                    X = IDAFT(X_QAM, c1, c2);
                    % X = IDAFT(X_QAM, 0, 0);
                    % 生成CPP
                    CPP = X(end - CPP_lenth + 1 : end) .* exp(-1j * 2 * pi * c1 * (N_data^2 + 2 * N_data * (-CPP_lenth:-1).'));
                    % CPP = X(end - CPP_lenth + 1 : end);
                    % 添加CPP / CP
                    S = [CPP; X];
                
                %% 叠加AWGN
                    Noise = sqrt(N0(i) / 2) * (randn(N_total,1) + 1j * randn(N_total,1));           % AWGN

                %% 接收端
                    % 接收到的时域符号
                    R = H * S + Noise;    % 接收信号 r = H * s + w
        
                    % 移除CPP
                    R = R(CPP_lenth + 1 : end);

                    % DAFT (c1 = c2 = 0时DAFT退化成为DFT)
                    Y = DAFT(R,c1,c2);
                    % Y = DAFT(R,0,0);

                    % 提取导频作信道估计
                    [H_eff,est_params] = Channel_Estimator(Y, pilot, pilot_idx, N_data, CPP_lenth, c1, c2, P);


                    % 均衡器
                    % Y_hat_full = Weighted_MRC_DFE(Y,H_eff,N0(i),active_idx,10);
                    Y_hat_full = MMSE(Y,H_eff(:,active_idx),N0(i));
                    Y_hat_active = Y_hat_full;

                    % 解映射
                    X0_hat = qamdemod(Y_hat_active, M, 'UnitAveragePower', true);
        
                    % 统计帧错误
                    [nErr,~] = biterr(X0_active,X0_hat);
                        
                    % 累加错误与统计总bit数
                    total_errs = total_errs + nErr;
                    total_bits = total_bits + N_active * bps;
                    % 实时显示进度，防止以为卡死
                    if mod(sym_count, 1e3) == 0
                        fprintf('已跑 %d 帧, 收集错误数: %d\n',sym_count, total_errs);
                    end
            end
            BER(i) = total_errs / total_bits;
            fprintf('SNR %d dB 完成: 共跑 %d 帧, BER = %.2e\n', SNRdB(i), sym_count, BER(i));
        end

%% 绘图结果展示
    figure(1);
    semilogy(SNRdB,BER,"-go","LineWidth",1.5);
    hold on;
    grid on;

%% 信道估计 - Demo
    % AFDM_Channel_Estimation_Debug_Demo();