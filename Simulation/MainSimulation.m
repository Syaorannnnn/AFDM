% 比较解耦后单导频 (EP) 与 CAZAC 序列导频在 MMSE 均衡下的 BER 性能

clear; clc; close all;

%% 基础参数与对象实例化

% --- 单导频 (EP) 配置 ---
configSingle = AfdmConfig();
configSingle.NumDataSubcarriers = 256;
configSingle.PilotType = "Single";
configSingle.PilotSnr = 35;

txSingle = AfdmTransmitter(configSingle);
rxSingle = AfdmReceiver(configSingle);
rxSingle.EqualizerType = "MMSE";
rxSingle.CsiMode = "Estimated";
rxSingle.EnableIterativeIc = false; % 单导频通常不使用迭代SIC

% --- CAZAC 导频配置 ---
configCazac = AfdmConfig();
configCazac.NumDataSubcarriers = 256;
configCazac.PilotType = "CAZAC";
configCazac.PilotSnr = 35;

txCazac = AfdmTransmitter(configCazac);
rxCazac = AfdmReceiver(configCazac);
rxCazac.EqualizerType = "MMSE";
rxCazac.CsiMode = "Estimated";
rxCazac.EnableIterativeIc = true; % CAZAC配合SIC消除干扰

%% 仿真环境设置
snrVectorDb = 0:5:30;
targetErrors = 1e2;
maxFrames = 1e5;

berSingle = zeros(length(snrVectorDb), 1);
berCazac = zeros(length(snrVectorDb), 1);

fprintf('AFDM BER 对比仿真开始\n');

%% 蒙特卡洛主循环
for snrIdx = 1:length(snrVectorDb)
    snrDb = snrVectorDb(snrIdx);

    noisePower = 10 ^ (-snrDb / 10);
    pilotPowerSingle = 10 ^ (configSingle.PilotSnr / 10) * noisePower;
    pilotPowerCazac = 10 ^ (configCazac.PilotSnr / 10) * noisePower;

    errSingle = 0; bitSingle = 0;
    errCazac = 0; bitCazac = 0;

    frameCount = 0;

    % 只要有一个方案还没跑够错误数，就继续跑（共享信道保证公平）
    while (errSingle < targetErrors || errCazac < targetErrors) && frameCount < maxFrames
        % 共享同一物理信道
        pathDelays = randperm(configSingle.MaxPathDelays + 1, configSingle.NumPaths) - 1;
        theta = (rand(1, configSingle.NumPaths) * 2 * pi) - pi;
        dopplers = configSingle.MaxNormDoppler * cos(theta);
        gains = (randn(1, configSingle.NumPaths) + 1j * randn(1, configSingle.NumPaths)) / sqrt(2 * configSingle.NumPaths);

        physicalChannelMatrix = LtvChannel(configSingle.TotalSubcarriers, pathDelays, dopplers, gains);

        % ==========================================
        % --- 链路 A: Single Pilot (EP) 处理 ---
        % ==========================================
        [txSigSingle, txDataSingle] = txSingle.transmit(pilotPowerSingle);
        % 生成独立噪声
        noiseSingle = sqrt(noisePower / 2) * (randn(size(txSigSingle)) + 1j * randn(size(txSigSingle)));
        rxSigSingle = physicalChannelMatrix * txSigSingle + noiseSingle;
        rxDataSingle = rxSingle.receive(rxSigSingle, noisePower, physicalChannelMatrix, pilotPowerSingle);

        [eS, ~] = biterr(txDataSingle, rxDataSingle);
        errSingle = errSingle + eS;
        bitSingle = bitSingle + configSingle.NumActiveCarriers * configSingle.BitsPerSymbol;

        % ==========================================
        % --- 链路 B: CAZAC Pilot 处理 ---
        % ==========================================
        [txSigCazac, txDataCazac] = txCazac.transmit(pilotPowerCazac);
        % 生成独立噪声
        noiseCazac = sqrt(noisePower / 2) * (randn(size(txSigCazac)) + 1j * randn(size(txSigCazac)));
        rxSigCazac = physicalChannelMatrix * txSigCazac + noiseCazac;
        rxDataCazac = rxCazac.receive(rxSigCazac, noisePower, physicalChannelMatrix, pilotPowerCazac);

        [eC, ~] = biterr(txDataCazac, rxDataCazac);
        errCazac = errCazac + eC;
        bitCazac = bitCazac + configCazac.NumActiveCarriers * configCazac.BitsPerSymbol;

        frameCount = frameCount + 1;

        if mod(frameCount, 1e2) == 0
            fprintf('CheckPoint : 已仿真 %4d 帧 | EP 错误数: %5d | CAZAC 错误数: %5d\r', frameCount, errSingle, errCazac);
        end

    end

    berSingle(snrIdx) = errSingle / bitSingle;
    berCazac(snrIdx) = errCazac / bitCazac;

    fprintf('SNR: %2d dB | 帧数: %4d | EP BER: %.4e | CAZAC BER: %.4e\n', ...
        snrDb, frameCount, berSingle(snrIdx), berCazac(snrIdx));
end

%% 绘图
figure;
semilogy(snrVectorDb, berSingle, '-bo', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy(snrVectorDb, berCazac, '-r^', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('AFDM BER: Single Pilot vs CAZAC (MMSE)');
legend('Single Pilot (EP) + MMSE', 'CAZAC Pilot + MMSE & SIC');
