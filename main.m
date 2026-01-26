clc; clear all; close all;
%% 基础参数
bps = 2;    % 比特率
M = 2^bps;  % 星座图点数

N_sym = 1e2;   % 总符号数
N_data = 2^8;  % 有效数据数
SNRdB = 0 : 5 : 20;  % 信噪比
SNR_idx_max = length(SNRdB);

P = 3; % 多径数
N0 = 10.^(-SNRdB ./ 10);
%% 多帧模拟
BER1 = zeros(1,SNR_idx_max);
BER2 = zeros(1,SNR_idx_max);
BER3 = zeros(1,SNR_idx_max);
BER4 = zeros(1,SNR_idx_max);
%% SNR步进
for i = 1 : SNR_idx_max
    % 初始化
    total_errs1 = 0;
    total_errs2 = 0;
    total_errs3 = 0;
    total_errs4 = 0;
    total_bits = N_data * N_sym * bps;
    for frame_idx = 1 : N_sym
    % maxNormDelay = 2;
    % maxNormDoppler = 1;

    % 给出固定的归一化时延、归一化多普勒、信道增益
    pathDelays = [0 1 2];
    pathDopplers = [-0.1 0.15 0.05];

    pathGains = (randn(1,P) + 1j .* randn(1,P)) ./ sqrt(2);   % 瑞利信道
    
    
    max_Doppler = max(abs(pathDopplers)); 
    max_Doppler_alpha = ceil(max_Doppler);              % 对多普勒向上取整
    max_Doppler_a = max_Doppler - max_Doppler_alpha;    % 多普勒的小数部分，范围[-0.5,0.5]
    max_Delay = max(pathDelays);
    % 均取前缀长度为最大多径时延
    CP_lenth = max_Delay;
    CPP_lenth = max_Delay;

    N = N_data + CPP_lenth; % 所需子载波个数 Data + Overhead
    data_range = (1 + CPP_lenth) : N;
    %% 设计c1、c2参数 
    k_v = 2;    % 用于对抗分数多普勒，为不同径带来足够高的分离度
    c1 = (2 * (max_Doppler_alpha + k_v) + 1) / (2 * N_data);           % 满足正交条件的c1，使多径各路径不发生交叠
    c2 = 1 / (N_data^2 * 2 * pi);
    
    if (2 * (max_Doppler + k_v) * (max_Delay + 1)) + max_Delay > N_data      %必须满足正交条件
        fprintf("子载波不满足正交！\n");
    end

        %% 发送端
        % 生成原始数据
        X0 = randi([0 M-1], N_data, 1);
        % 映射成星座点
        X_QAM = qammod(X0, M, 'UnitAveragePower', true);
        % IDAFT(c1 = c2 = 0时IDAFT退化成为IDFT)
        X = IDAFT(X_QAM, c1, c2);
        X_OFDM = IDAFT(X_QAM,0,0); 
        % 生成CPP
        CPP = X(end - CPP_lenth + 1 : end) .* exp(-1j * 2 * pi * c1 * (N_data^2 + 2 * N_data * (-CPP_lenth:-1).'));
        CP = X_OFDM(end - CP_lenth + 1 : end);
        % CPP = X(end - CPP_lenth + 1 : end);
        % 添加CPP / CP
        S = [CPP; X];
        S_OFDM = [CP; X_OFDM];
        
        %% 经过LTV信道叠加AWGN噪声
        Noise = sqrt(N0(i) / 2) * (randn(N,1) + 1j * randn(N,1));           % 构造AWGN噪声
        H = LTVchannel(N,pathDelays,pathDopplers,pathGains);                % 构造LTV信道矩阵
        H_eff = Gen_Eff_Channel(H,CPP_lenth,c1,c2);                         % AFDM的有效信道矩阵
        H_eff_OFDM = Gen_Eff_Channel(H,CP_lenth,0,0);                       % OFDM的有效信道矩阵

        %% 接收端
        R_AFDM = H * S + Noise;    % 接收信号 r = H * s + w
        R_OFDM = H * S_OFDM + Noise;
        % 移除CPP / CP
        R_AFDM = R_AFDM(CPP_lenth + 1 : end);
        R_OFDM = R_OFDM(CP_lenth + 1 : end);
        % DAFT (c1 = c2 = 0时DAFT退化成为DFT)
        Y_AFDM = DAFT(R_AFDM,c1,c2);
        Y_OFDM = DAFT(R_OFDM,0,0);
        % MMSE / MRC-Based DFE
        Y_AFDM_hat1 = MMSE(Y_AFDM,H_eff,N0(i));
        Y_AFDM_hat2 = Weighted_MRC_DFE(Y_AFDM,H_eff,N0(i),6);

        Y_OFDM_hat1 = MMSE(Y_OFDM,H_eff_OFDM,N0(i));
        Y_OFDM_hat2 = Weighted_MRC_DFE(Y_OFDM,H_eff_OFDM,N0(i),6);
        % 解映射
        X0_hat1 = qamdemod(Y_AFDM_hat1, M, 'UnitAveragePower', true);
        X0_hat2 = qamdemod(Y_AFDM_hat2, M, 'UnitAveragePower', true);
        X0_hat3 = qamdemod(Y_OFDM_hat1, M, 'UnitAveragePower', true);
        X0_hat4 = qamdemod(Y_OFDM_hat2, M, 'UnitAveragePower', true);

        % % 1. 检查发送和接收能量是否一致 (排查归一化问题)
        % disp(['Tx 能量: ', num2str(norm(S_OFDM)), ' | Rx 能量: ', num2str(norm(R_OFDM))]);
        % 
        % % 2. 检查 H_eff 是否真的是单位阵 (排查信道矩阵生成问题)
        % % 理想情况下，H_eff_OFDM 应该对角线全是 1，其余全为 0
        % disp(['H_eff 对角线均值: ', num2str(mean(abs(diag(H_eff_OFDM))))]);
        % disp(['H_eff 非对角线最大值: ', num2str(max(max(abs(H_eff_OFDM - diag(diag(H_eff_OFDM))))))]);
        % 
        % % 3. 直接对比频域信号 (核心判决！)
        % % 如果 Y_OFDM 和 X_QAM 不相等，说明 DAFT/IDAFT 不匹配或时域截取错位
        % diff_Y = norm(Y_OFDM - X_QAM);
        % disp(['Y_OFDM 与 X_QAM 的误差: ', num2str(diff_Y)]);
        % 
        % % 4. 直接对比均衡后的信号
        % diff_Xhat = norm(X_OFDM_hat1 - X_QAM);
        % disp(['均衡后与原信号的误差: ', num2str(diff_Xhat)]);

        % 统计帧错误
        [nErr1,~] = biterr(X0,X0_hat1);
        [nErr2,~] = biterr(X0,X0_hat2);
        [nErr3,~] = biterr(X0,X0_hat3);
        [nErr4,~] = biterr(X0,X0_hat4);
        
        % 累加错误与总bit数
        total_errs1 = total_errs1 + nErr1;
        total_errs2 = total_errs2 + nErr2;
        total_errs3 = total_errs3 + nErr3;
        total_errs4 = total_errs4 + nErr4;

    end
    BER1(i) = total_errs1 / total_bits;
    BER2(i) = total_errs2 / total_bits;
    BER3(i) = total_errs3 / total_bits;
    BER4(i) = total_errs4 / total_bits;
    BER1(BER1 == 0) = 1e-6;
    BER2(BER2 == 0) = 1e-6;
    BER3(BER3 == 0) = 1e-6;
    BER4(BER4 == 0) = 1e-6;
end
%% 
% fprintf("BER via MMSE %.3e\nBER via MRC_DFE %.3e\n",BER1,BER2);
set(0, 'defaultfigurecolor', 'w');
Amp1 = abs(H_eff);
Amp2 = abs(H_eff_OFDM);
Amp = abs(H);
figure(1);
set(gcf,'position',[250 300 600 200]);
subplot(1,3,1);
imagesc(Amp1);
set(gca,'position',[0.05 0.3 0.3 0.6]);
axis xy;
colorbar;
title('AFDM输入输出等效矩阵Heff示意图');
xlabel('发送信号子载波索引');
ylabel('接收信号子载波索引')
subplot(1,3,2);
imagesc(Amp2);
set(gca,'position',[0.36 0.3 0.3 0.6]);
axis xy;
colorbar;
title('OFDM输入输出等效矩阵Heff示意图');
xlabel('发送信号子载波索引');
ylabel('接收信号子载波索引');
subplot(1,3,3);
imagesc(Amp);
set(gca,'position',[0.7 0.3 0.3 0.6]);
axis xy;
colorbar;
title('LTV信道矩阵示意图');
xlabel('发送信号子载波索引');
ylabel('接收信号子载波索引');


figure(2);
semilogy(SNRdB,BER1,"-go",SNRdB,BER2,"--b^",SNRdB,BER3,"-cs",SNRdB,BER4,"--rd","LineWidth",1.5);
% semilogy(SNRdB,BER1,"-go",SNRdB,BER2,"--b^","LineWidth",1.5);
title("LMMSE vs MRC-Based DFE","BER");
ylim([1e-6,1]);
xlabel("SNR in dB");
ylabel("BER");
legend("AFDM,MMSE","AFDM,MRC-Based DFE","OFDM,MMSE","OFDM,MRC-Based DFE");
hold on;
grid on;

%% M-QAM理论BER
% BER = zeros(1,SNR_idx_max);
% for i = 1 : SNR_idx_max
%     X = randi([0 M-1], 1e3, 1);
%     X_111 = qammod(X, M, 'UnitAveragePower', true);
%     Y = qamdemod(awgn(X_111,SNRdB(i)),M, 'UnitAveragePower', true);
%     [~,BER(i)] = biterr(X,Y);
%     BER(BER == 0) = 1e-6;
% end