%
%   功能： 仿真主测试界面
%
clc; 
clear; 
close all;
%% 仿真参数
carrierFreq = 4e9;  % 载波频率
bandwidth = 1e6;    % 带宽
samplingFreq = 2 * bandwidth; % 奈奎斯特采样频率

bps = 2;    % 比特率
M = 2^bps;  % 星座图点数
N = 2^8;    % 子载波个数
NT = 1;     % 符号周期数

%% Channel Detection 多帧信号经过DD信道叠加高斯白噪
N_frames = 100;                         % 总帧数
SNRdB = 0 : 3 : 30 ;                    % 预设最大信噪比30dB
SNR_idx_max = length(SNRdB);
% 构建双色散信道
P = 3;
max_normDelay = 2;           % 最大归一化延迟索引
max_normDoppler = 1;         % 最大归一化多普勒

BER1 = zeros(1,SNR_idx_max);
BER2 = zeros(1,SNR_idx_max);
for SNR_idx = 1 : SNR_idx_max
    N0 = 10^(-SNRdB(SNR_idx) / 10);                  % 噪声功率
    % Window_Power_Scale = 1 / 0.375;
    % N0 = N0 * Window_Power_Scale;

    % 初始化误差统计
    total_errors_LMMSE = 0;
    total_errors_MRC_DFE = 0;
    total_valid_bits = 0;

        for frame_idx = 1 : N_frames
            % pathDelays = randi(max_normDelay, [1,P]);
            % pathDelays = sort(pathDelays - min(pathDelays));      % 整数
            % pathDopplers = max_normDoppler - 2 * max_normDoppler * rand(1,P);

            pathDelays = [0 1 2];
            pathDopplers = [-0.8 0.5 0.4];

            pathGains = (randn(1,P) + 1j .* randn(1,P)) ./ sqrt(2);   % 瑞利信道

            max_Doppler = max(abs(pathDopplers));
            max_Doppler_alpha = ceil(max_Doppler);              % 对多普勒向上取整
            max_Doppler_a = max_Doppler - max_Doppler_alpha;    % 多普勒的小数部分，范围[-0.5,0.5]
            max_Delay = max(pathDelays);

            % 设计c1、c2参数 
            k_v = 2;    % 用于对抗分数多普勒带来的IDI与IDPI，为不同径带来足够高的分离度
            c1_AFDM = (2 * (max_Doppler_alpha + k_v) + 1) / (2 * N);           % 满足正交条件的c1，使多径各路径不发生交叠
            c2_AFDM = 1 / (N^2 * 2 * pi);

            if (2 * (max_Doppler + k_v) * (max_Delay + 1)) + max_Delay > N      %必须满足正交条件
                fprintf("子载波不满足正交！\n");
            end
        
            guardWidth = (max_Delay + 1) * (2 * (max_Doppler_alpha + k_v) + 1) - 1;     % 保护间隔宽度
            Pilot_lenth = max_Delay + 1;                                                  % Pilot_lenth >= maxDelay - 1
            Data_lenth = N - 2 * guardWidth - Pilot_lenth;
        
            H = DDchannel(N,pathDelays,pathDopplers,pathGains,samplingFreq,c1_AFDM);
            H_eff = AFDM_Channel(H,c1_AFDM,c2_AFDM);

            X_i = randi([0 M-1], Data_lenth, NT);
            % 生成能量归一化的M-QAM复符号
            X = qammod(X_i, M, 'UnitAveragePower',true); 
            
            % 构造AFDM帧
            % 构造独立的随机 Pilot
            Pilot_i = randi([0 M-1], Pilot_lenth, 1);
            Pilot = qammod(Pilot_i, M, 'UnitAveragePower', true);     % Pilot
            % 帧结构：[ Pilot | Null_Guard | Data | Null_Guard ]
            X_AFDM = [Pilot;zeros(guardWidth,1);X;zeros(guardWidth,1)];
            X_Pilot_Only = [Pilot;zeros(N - Pilot_lenth,1)];
            % 有用信号范围 
            data_range = (Pilot_lenth + guardWidth + 1) : (Pilot_lenth + guardWidth + Data_lenth);
            % theta = Embbed_Pilot_Estimation(X_AFDM,N,c1_AFDM,c2_AFDM,P,max_Delay,max_Doppler);

            % 截短矩阵T
            I_N = eye(N);
            % T = I_N((guardWidth - (max_Doppler_alpha + k_v)) : N - (max_Doppler_alpha + k_v) - 1, :);
            T = I_N(data_range, :);

            Noise = sqrt(N0 / 2) * (randn(N,1) + 1j * randn(N,1));  % 构造AWGN噪声
            Y_AFDM = H_eff * X_AFDM + AFDM_Demodulate(Noise,c1_AFDM,c2_AFDM);    %接收信号 y = H_eff * x + Aw
            Y_PI = H_eff * X_Pilot_Only;
            Y_Detec = Y_AFDM - Y_PI;
            H_eff_truncated = H_eff * T';
            % search_count = P * 4; 
            % pathRows = Find_PathRows(H_eff_truncated, search_count);
            % X_truncated = T * X_AFDM;
            % 比较LMMSE MRC-DFE
            [X_hat1] = Iterative_MMSE(Y_AFDM,H_eff_truncated,N0,3);
            [X_hat2] = Weighted_MRC_DFE(Y_AFDM,H_eff,data_range,N0,3);

            RX_X1 = qamdemod(X_hat1, M, 'UnitAveragePower', true);
            RX_X2 = qamdemod(X_hat2, M, 'UnitAveragePower', true);

            %统计本帧错误
            [nErr1, ~] = biterr(X_i,RX_X1);
            [nErr2, ~] = biterr(X_i,RX_X2);

            total_errors_LMMSE = total_errors_LMMSE + nErr1;
            total_errors_MRC_DFE = total_errors_MRC_DFE + nErr2;
            total_valid_bits = total_valid_bits + Data_lenth;

        end

        BER1(SNR_idx) = total_errors_LMMSE / total_valid_bits;
        if(BER1(SNR_idx) == 0)
            BER1(SNR_idx) = 1e-15;
        end
        BER2(SNR_idx) = total_errors_MRC_DFE / total_valid_bits;
        if(BER2(SNR_idx) == 0)
            BER2(SNR_idx) = 1e-15;
        end

end


%% 绘图
% Amp = abs(H_eff);
% Amp1 = abs(H_Windowed);
% figure(1);
% subplot(1,2,1)
% imagesc(Amp);
% colorbar;
% colormap(slanCM('rainbow'));
% title('AFDM输入输出矩阵Heff示意图')
% xlabel('发送信号子载波索引')
% ylabel('接收信号子载波索引')
% subplot(1,2,2)
% imagesc(Amp1);
% colorbar;
% title('重构Heff示意图')
% xlabel('发送信号子载波索引')
% ylabel('接收信号子载波索引')
% 
% figure(2);
% [X1, Y1] = meshgrid(1 : size(Amp, 1), 1 : size(Amp, 2)); % 生成网格点
% surf(X1,Y1,Amp); % 绘制曲面图
% colormap(slanCM('rainbow'))
% shading interp;
% colorbar;
% title('AFDM输入输出矩阵H_e_f_f示意图')
% xlabel('接收信号子载波索引')
% ylabel('发送信号子载波索引')
% zlabel('幅度')

% figure(3);
% stem(abs(Y_AFDM),'o');
% title('接收信号幅度示意图')
% xlabel('接收信号子载波索引')
% ylabel('幅值')

% figure(4);
% subplot(1,2,1);

% plot(1:1:N - guardWidth,x1);
% xlim([0,N - 1]);
% ylim([0,15]);
% subplot(1,2,2);

% plot(1:1:N - guardWidth,x3);
% xlim([0,N - 1]);
% ylim([0,15]);

% scatterplot(X);
% scatterplot(X_hat1);
% title('LMMSE');
% scatterplot(X_hat2);
% title('MRC-DFE');

% fprintf('SER via LMMSE: %.4e \nSER via MRC-DFE: %.4e \n',ser1,ser2);
%%
figure;
semilogy(SNRdB,BER1,"--go",SNRdB,BER2,"--b^","LineWidth",1.5);
title("LMMSE与MRC-Based DFE BER比较");
xlabel("SNR in dB");
ylabel("BER");
legend("LMMSE","MRC-Based DFE");
hold on;
grid on;