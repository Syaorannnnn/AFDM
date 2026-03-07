% 比较解耦后单导频 (EP) 与 CAZAC 序列导频在 MMSE 均衡下的 BER 性能

clear; clc; close all;

%% 启动并行池
if isempty(gcp('nocreate'))
    parpool(5); 
end

%% 全局基础参数定义

numDataSubcarriers = 256;
pilotSnr = 35; 

% 建立外部参考对象仅用于获取衍生参数(如有效子载波数量)，坚决不传进 parfor
configSingleRef = AfdmConfig();
configSingleRef.NumDataSubcarriers = numDataSubcarriers;
configSingleRef.PilotType = "Single";

configCazacRef = AfdmConfig();
configCazacRef.NumDataSubcarriers = numDataSubcarriers;
configCazacRef.PilotType = "CAZAC";

% [提取纯量参数] 供 parfor 内部信道生成和比特数统计使用
maxPathDelays = configSingleRef.MaxPathDelays;
numPaths = configSingleRef.NumPaths;
maxNormDoppler = configSingleRef.MaxNormDoppler;
totalSubcarriers = configSingleRef.TotalSubcarriers;

numActiveCarriersSingle = configSingleRef.NumActiveCarriers;
numActiveCarriersCazac = configCazacRef.NumActiveCarriers;
bitsPerSymbol = configSingleRef.BitsPerSymbol;

%% 仿真环境设置
snrVectorDb = 0:5:30;
targetErrors = 1e3;
maxFrames = 1e6;
batchSize = 1e2; 

berSingle = zeros(length(snrVectorDb), 1);
berCazac = zeros(length(snrVectorDb), 1);

fprintf('AFDM BER 对比仿真开始\n');

%% 蒙特卡洛主循环
for snrIdx = 1:length(snrVectorDb)
    snrDb = snrVectorDb(snrIdx);

    noisePower = 10 ^ (-snrDb / 10);
    pilotPower = 10 ^ (pilotSnr / 10) * noisePower;

    errSingle = 0; bitSingle = 0;
    errCazac = 0; bitCazac = 0;

    frameCount = 0;

    while (errSingle < targetErrors || errCazac < targetErrors) && frameCount < maxFrames
        
        currentBatchSize = min(batchSize, maxFrames - frameCount);
        
        batchErrSingle = zeros(1, currentBatchSize);
        batchErrCazac  = zeros(1, currentBatchSize);

        parfor i = 1:currentBatchSize
            % ==========================================================
            % [核心优化]：在 Worker 内部按需实例化 Handle 对象
            % MATLAB 对象创建开销极低，此举彻底解决了广播警告与内部状态修改的数据竞争
            % ==========================================================
            
            % --- 单导频 (EP) 链路局部实例化 ---
            localConfigSingle = AfdmConfig();
            localConfigSingle.NumDataSubcarriers = numDataSubcarriers;
            localConfigSingle.PilotType = "Single";
            localConfigSingle.PilotSnr = pilotSnr;
            
            localTxSingle = AfdmTransmitter(localConfigSingle);
            localRxSingle = AfdmReceiver(localConfigSingle);
            localRxSingle.EqualizerType = "MMSE";
            localRxSingle.CsiMode = "Estimated";
            localRxSingle.EnableIterativeIc = false;
            
            % --- CAZAC 导频链路局部实例化 ---
            localConfigCazac = AfdmConfig();
            localConfigCazac.NumDataSubcarriers = numDataSubcarriers;
            localConfigCazac.PilotType = "CAZAC";
            localConfigCazac.PilotSnr = pilotSnr;
            
            localTxCazac = AfdmTransmitter(localConfigCazac);
            localRxCazac = AfdmReceiver(localConfigCazac);
            localRxCazac.EqualizerType = "MMSE";
            localRxCazac.CsiMode = "Estimated";
            localRxCazac.EnableIterativeIc = true;

            % --- 独立随机信道生成 ---
            pathDelays = randperm(maxPathDelays + 1, numPaths) - 1;
            theta = (rand(1, numPaths) * 2 * pi) - pi;
            dopplers = maxNormDoppler * cos(theta);
            gains = (randn(1, numPaths) + 1j * randn(1, numPaths)) / sqrt(2 * numPaths);

            physicalChannelMatrix = LtvChannel(totalSubcarriers, pathDelays, dopplers, gains);

            % --- 链路 A: Single Pilot (EP) 发送与接收 ---
            [txSigSingle, txDataSingle] = localTxSingle.transmit(pilotPower);
            noiseSingle = sqrt(noisePower / 2) * (randn(size(txSigSingle)) + 1j * randn(size(txSigSingle)));
            rxSigSingle = physicalChannelMatrix * txSigSingle + noiseSingle;
            rxDataSingle = localRxSingle.receive(rxSigSingle, noisePower, physicalChannelMatrix, pilotPower);

            [batchErrSingle(i), ~] = biterr(txDataSingle, rxDataSingle);

            % --- 链路 B: CAZAC Pilot 发送与接收 ---
            [txSigCazac, txDataCazac] = localTxCazac.transmit(pilotPower);
            noiseCazac = sqrt(noisePower / 2) * (randn(size(txSigCazac)) + 1j * randn(size(txSigCazac)));
            rxSigCazac = physicalChannelMatrix * txSigCazac + noiseCazac;
            rxDataCazac = localRxCazac.receive(rxSigCazac, noisePower, physicalChannelMatrix, pilotPower);

            [batchErrCazac(i), ~] = biterr(txDataCazac, rxDataCazac);

            delete(localTxSingle);
            delete(localRxSingle);
            delete(localConfigSingle);
            
            delete(localTxCazac);
            delete(localRxCazac);
            delete(localConfigCazac);
        end

        % 数据归集
        errSingle = errSingle + sum(batchErrSingle);
        errCazac  = errCazac + sum(batchErrCazac);

        frameCount = frameCount + currentBatchSize;
        
        bitSingle = bitSingle + currentBatchSize * numActiveCarriersSingle * bitsPerSymbol;
        bitCazac  = bitCazac + currentBatchSize * numActiveCarriersCazac * bitsPerSymbol;

        fprintf('CheckPoint : 已仿真 %4d 帧 | EP 错误数: %5d | CAZAC 错误数: %5d\r', frameCount, errSingle, errCazac);
    end
    
    fprintf('\n'); 

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
