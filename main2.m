clc; clear all; close all;
%% 基础参数
bps = 2;    % 比特率
M = 2^bps;  % 星座图点数
N_data = 2^8;  % 有效数据数
SNRdB = 0 : 5 : 20;  % 信噪比
SNR_idx_max = length(SNRdB);

P = 3; % 多径数
N0 = 10.^(-SNRdB ./ 10);

target_errors = 100;
max_sym = 1e6;

maxNormDoppler = 2;
% 给出固定的归一化时延，初始化归一化多普勒、信道增益
pathDelays = [0 1 2];
pathDopplers = zeros(P);
pathGains = zeros(P);

max_Delay = max(pathDelays);

CPP_lenth = max_Delay;

N = N_data + CPP_lenth; % 所需子载波个数 Data + Overhead
data_range = (1 + CPP_lenth) : N;

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
            while(total_errs <= target_errors && sym_count < max_sym)
                % Jakes分布
                % P个随机相位theta
                theta = (rand(1, P) * 2 * pi) - pi;
                pathDopplers = maxNormDoppler * cos(theta);
                
                % 瑞利信道
                pathGains = (randn(1, P) + 1j * randn(1, P)) / sqrt(2 * P);

                max_Doppler = max(abs(pathDopplers)); 
                max_Doppler_alpha = round(max_Doppler);              % 对多普勒向上取整
                max_Doppler_a = max_Doppler - max_Doppler_alpha;    % 多普勒的小数部分，范围[-0.5,0.5]
                %% 设计c1、c2参数 
                k_v = 2;    % 用于对抗分数多普勒，为不同径带来足够高的分离度
                c1 = (2 * (max_Doppler_alpha + k_v) + 1) / (2 * N_data);           % 满足正交条件的c1，使多径各路径不发生交叠
                c2 = 1 / (N_data^2 * 2 * pi);
                
                if (2 * (max_Doppler + k_v) * (max_Delay + 1)) + max_Delay > N_data      %必须满足正交条件
                    fprintf("子载波不满足正交！\n");
                end

                Q_calc = (max_Delay + 1) * (2 * (maxNormDoppler + k_v) + 1) - 1;
                Q = ceil(Q_calc);
                if Q >= N_data
                    fprintf("ZP长度过大,没有资源发数据了!请增加N_data或减小多普勒/时延\n");
                end
                
                N_active = N_data - Q; % 有效数据子载波数
                active_idx = 1: N_active;   % 为方便起见，先把ZP全部加在有效数据末尾

                H = LTVchannel(N,pathDelays,pathDopplers,pathGains);                % 构造LTV信道矩阵
                H_eff = Gen_Eff_Channel(H,CPP_lenth,c1,c2);                         % AFDM的有效信道矩阵
                % H_eff = Gen_Eff_Channel(H,CPP_lenth,0,0);                         % OFDM的有效信道矩阵

                % 每生成一个符号，计数器+1
                sym_count = sym_count + 1;
    
                %% 发送端
                % 生成原始数据
                X0_active = randi([0 M-1], N_active, 1);
                % 映射成星座点
                X_QAM_active = qammod(X0_active, M, 'UnitAveragePower', true);
                % 组建完整的DAFT域符号向量（加入ZP）
                X_QAM = zeros(N_data,1);
                X_QAM(active_idx) = X_QAM_active;
                % IDAFT(c1 = c2 = 0时IDAFT退化成为IDFT)
                X = IDAFT(X_QAM, c1, c2);
                % X = IDAFT(X_QAM, 0, 0);
                % 生成CPP
                CPP = X(end - CPP_lenth + 1 : end) .* exp(-1j * 2 * pi * c1 * (N_data^2 + 2 * N_data * (-CPP_lenth:-1).'));
                % CPP = X(end - CPP_lenth + 1 : end);
                % 添加CPP / CP
                S = [CPP; X];
                
                %% 经过LTV信道叠加AWGN
                Noise = sqrt(N0(i) / 2) * (randn(N,1) + 1j * randn(N,1));           % AWGN
    
                %% 接收端
                R = H * S + Noise;    % 接收信号 r = H * s + w
    
                % 移除CPP
                R = R(CPP_lenth + 1 : end);

                % DAFT (c1 = c2 = 0时DAFT退化成为DFT)
                Y = DAFT(R,c1,c2);
                % Y = DAFT(R,0,0);

                % MMSE / MRC-Based DFE
                Y_hat_full = Weighted_MRC_DFE(Y,H_eff,N0(i),active_idx,5);
                % Y_hat_full = MMSE(Y,H_eff(:,active_idx),N0(i));
                Y_hat_active = Y_hat_full(active_idx);
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

figure(1);
semilogy(SNRdB,BER,"-go","LineWidth",1.5);
hold on;
grid on;