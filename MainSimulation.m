function MainSimulation()
    OfdmVersusAfdmStimulation();
end

% --- OFDM vs AFDM ---
function OfdmVersusAfdmStimulation()
    clc; clear; close all;
    % --- 仿真全局参数设置 ---
    snrRange = 0 : 5 : 20;           
    % 建议将 targetErrors 调大至 1e4，以确保高信噪比下的曲线平滑度
    targetErrors = 1e3;           
    maxFrames = 1e7;            

    waveforms = ["OFDM", "AFDM"];
    equalizers = ["MMSE", "MRC-DFE"];

    % 定义线型和颜色
    lineStyles = {"--b", "--r"; "-g", "-m"}; 
    markers = {"o", "s"; "^", "d"};

    % 预分配结果存储矩阵 (维度: 波形 x 均衡器 x SNR)
    berResults = zeros(length(waveforms), length(equalizers), length(snrRange));

    % 实例化系统
    sys = AfdmSystem();
    sys.CsiMode = "Perfect";
    sys.NumMaxIterations = 8;

    % --- 绘图准备 ---
    figure('Name', 'Waveform & Equalizer Comparison', 'Position', [100, 100, 800, 600]);
    hold on; grid on;
    set(gca, 'YScale', 'log');
    xlabel('SNR_d(dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    title('AFDM vs OFDM under Doubly Dispersive Channel', 'FontSize', 14);

    % --- 主循环开始 ---
    for sIdx = 1:length(snrRange)
        currentSnr = snrRange(sIdx);
        noisePower = 10^(-currentSnr / 10);
        % pilotPower = 10^(sys.PilotSnr / 10) * noisePower;
        
        fprintf('========================================\n');
        fprintf('正在运行 SNR: %d dB...\n', currentSnr);
        
        % 初始化当前 SNR 下的错误计数器和总比特数
        totalErrors = zeros(length(waveforms), length(equalizers));
        totalBits = zeros(length(waveforms), 1);
        frameCount = 0;
        
        % 所有组合的错误数都达到targetErrors再停止统计, 以确保最高精度的统计
        while min(totalErrors(:)) < targetErrors && frameCount < maxFrames
            
            % 生成全局共享的物理信道参数 (控制变量法的核心)
            % 确保同一个 Frame 下，OFDM 和 AFDM 经历相同的衰落特征
            pathDelays = randperm(sys.MaxPathDelays + 1, sys.NumPaths) - 1;
            theta = (rand(1, sys.NumPaths) * 2 * pi) - pi;
            dopplers = sys.MaxNormDoppler * cos(theta);
            gains = (randn(1, sys.NumPaths) + 1j * randn(1, sys.NumPaths)) / sqrt(2 * sys.NumPaths);
            
            % % 发送前修正导频功率
            % sys.PilotSymbol = sqrt(pilotPower) * (sys.PilotSymbol / abs(sys.PilotSymbol));
            
            % 2. 遍历不同波形
            for wIdx = 1:length(waveforms)
                sys.WaveformType = waveforms(wIdx);
                
                % 使用当前帧的共享参数生成物理信道矩阵
                physicalChannelMatrix = LtvChannel(sys.TotalSubcarriers, pathDelays, dopplers, gains);
                
                % 发送端生成信号
                [txSignal, txData] = sys.transmit();
                
                % 通过信道并加入噪声 (为当前波形生成一次独立的噪声)
                noise = sqrt(noisePower/2) * (randn(size(txSignal)) + 1j * randn(size(txSignal)));
                rxSignal = physicalChannelMatrix * txSignal + noise;
                
                % 3. 遍历不同均衡器
                % 关键点：MMSE 和 MRC-DFE 此时处理的是完全相同的 rxSignal 和 physicalChannelMatrix
                for eIdx = 1:length(equalizers)
                    sys.EqualizerType = equalizers(eIdx);
                    
                    % 接收与均衡
                    rxData = sys.receive(rxSignal, noisePower, physicalChannelMatrix);
                    
                    % 统计当前组合的误码率
                    [errs, ~] = biterr(txData, rxData);
                    totalErrors(wIdx, eIdx) = totalErrors(wIdx, eIdx) + errs;
                end
                
                % 累计该波形下的总传输比特数
                totalBits(wIdx) = totalBits(wIdx) + sys.NumActiveCarriers * sys.BitsPerSymbol;
            end
            
            frameCount = frameCount + 1;
        end
        
        % 计算并记录当前 SNR 的最终 BER
        for wIdx = 1:length(waveforms)
            for eIdx = 1:length(equalizers)
                berResults(wIdx, eIdx, sIdx) = totalErrors(wIdx, eIdx) / totalBits(wIdx);
            end
        end
        
        % 打印当前 SNR 进度
        fprintf('  -> 完成. 帧数: %d\n', frameCount);
        fprintf('     OFDM  MMSE: %.2e | MRC-DFE: %.2e\n', berResults(1, 1, sIdx), berResults(1, 2, sIdx));
        fprintf('     AFDM  MMSE: %.2e | MRC-DFE: %.2e\n', berResults(2, 1, sIdx), berResults(2, 2, sIdx));
    end

    % --- 统一绘制最终曲线 ---
    legendEntries = strings(1, length(waveforms) * length(equalizers));
    plotIndex = 1;

    for wIdx = 1:length(waveforms)
        for eIdx = 1:length(equalizers)
            semilogy(snrRange, squeeze(berResults(wIdx, eIdx, :)), ...
                    lineStyles{wIdx, eIdx}, 'Marker', markers{wIdx, eIdx}, ...
                    'LineWidth', 2, 'MarkerSize', 8);
            legendEntries(plotIndex) = sprintf('%s - %s', waveforms(wIdx), equalizers(eIdx));
            plotIndex = plotIndex + 1;
        end
    end

    legend(legendEntries, 'Location', 'southwest', 'FontSize', 11);
    fprintf('全部仿真完成！\n');
end