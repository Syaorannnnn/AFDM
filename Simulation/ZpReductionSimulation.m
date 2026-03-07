% ZpReductionSimulation.m
% 核心实验: 缩减 ZP 下 Single EP 与 CAZAC 导频的 BER 对比
%
% 实验设计:
%   CAZAC 序列导频的核心优势在于提供更鲁棒的信道估计,
%   使系统可以在缩减保护带 (ZP) 的条件下仍保持可接受的检测性能,
%   从而提升频谱效率。
%
%   本仿真将 ZP 从理论完整值逐步缩减, 观察两种导频策略的 BER 衰退速度。
%   预期: Single EP 在 ZP 缩减时 BER 快速恶化,
%         CAZAC + SIC 的 BER 衰退更平缓。
%
% 依赖文件:
%   AfdmConfig, AfdmTransmitter, AfdmReceiver, AfdmTransforms, LtvChannel

clear; clc; close all;

%% ====================================================================
%  第一部分: 仿真参数
%  ====================================================================

% ZP 扫描: 从理论完整值 (77) 逐步缩减
zpValues = [77, 55, 40, 30, 20, 15, 10];

% 选取有代表性的 SNR 点
snrValues = [10, 15, 20];

% 蒙特卡洛参数
targetErrors = 5e2; % 每个 SNR 点至少累积的错误比特数
maxFrames = 1e5; % 单个 SNR 点的最大帧数上限

%% ====================================================================
%  第二部分: 预览实验配置
%  ====================================================================

numZp = length(zpValues);
numSnr = length(snrValues);

berSingle = zeros(numZp, numSnr);
berCazac = zeros(numZp, numSnr);
overheadSingle = zeros(numZp, 1);
overheadCazac = zeros(numZp, 1);
activeCarriersSingle = zeros(numZp, 1);
activeCarriersCazac = zeros(numZp, 1);

fprintf('========================================================\n');
fprintf('  CAZAC vs Single EP: 缩减 ZP 下的 BER 对比实验\n');
fprintf('========================================================\n');
fprintf('ZP 扫描值: '); fprintf('%d  ', zpValues);
fprintf('\nSNR 扫描值: '); fprintf('%d  ', snrValues); fprintf('dB\n');
fprintf('目标错误数: %d, 最大帧数: %d\n\n', targetErrors, maxFrames);

fprintf('%6s | %12s | %12s | %12s | %12s\n', 'ZP', 'EP开销', 'EP有效载波', 'CAZAC开销', 'CAZAC有效载波');

for zpIdx = 1:numZp
    zp = zpValues(zpIdx);
    spanEP = 2 * zp + 1; % Single: 导频长度=1
    spanCA = 2 * zp + 15; % CAZAC: 导频长度=15
    activeEP = 256 - spanEP;
    activeCA = 256 - spanCA;
    fprintf('%6d | %10.1f%% | %12d | %10.1f%% | %12d\n', ...
        zp, spanEP / 256 * 100, activeEP, spanCA / 256 * 100, activeCA);
end

fprintf('\n');

%% ====================================================================
%  第三部分: 蒙特卡洛主循环
%  ====================================================================

for zpIdx = 1:numZp
    currentZp = zpValues(zpIdx);

    % --- Single EP 配置 ---
    cfgSingle = AfdmConfig();
    cfgSingle.PilotType = "Single";
    cfgSingle.ManualZeroPaddingLength = currentZp;

    txSingle = AfdmTransmitter(cfgSingle);
    rxSingle = AfdmReceiver(cfgSingle);
    rxSingle.EqualizerType = "MMSE";
    rxSingle.CsiMode = "Estimated";
    rxSingle.EnableIterativeIc = false;

    % --- CAZAC 配置 ---
    cfgCazac = AfdmConfig();
    cfgCazac.PilotType = "CAZAC";
    cfgCazac.ManualZeroPaddingLength = currentZp;

    txCazac = AfdmTransmitter(cfgCazac);
    rxCazac = AfdmReceiver(cfgCazac);
    rxCazac.EqualizerType = "MMSE";
    rxCazac.CsiMode = "Estimated";
    rxCazac.EnableIterativeIc = true;
    % SIC 使用自适应收敛 (默认: 最大 8 轮, 阈值 1e-3)

    % 记录开销
    overheadSingle(zpIdx) = 1 - cfgSingle.NumActiveCarriers / cfgSingle.NumDataSubcarriers;
    overheadCazac(zpIdx) = 1 - cfgCazac.NumActiveCarriers / cfgCazac.NumDataSubcarriers;
    activeCarriersSingle(zpIdx) = cfgSingle.NumActiveCarriers;
    activeCarriersCazac(zpIdx) = cfgCazac.NumActiveCarriers;

    for snrIdx = 1:numSnr
        snrDb = snrValues(snrIdx);
        noisePower = 10 ^ (-snrDb / 10);
        pilotPowerSingle = 10 ^ (cfgSingle.PilotSnr / 10) * noisePower;
        pilotPowerCazac = 10 ^ (cfgCazac.PilotSnr / 10) * noisePower;

        errS = 0; bitS = 0;
        errC = 0; bitC = 0;
        frameCount = 0;

        while (errS < targetErrors || errC < targetErrors) && frameCount < maxFrames

            % 共享物理信道 (分数多普勒)
            pathDelays = randperm(cfgSingle.MaxPathDelays + 1, cfgSingle.NumPaths) - 1;
            dopplers = cfgSingle.MaxNormDoppler * (2 * rand(1, cfgSingle.NumPaths) - 1);
            gains = (randn(1, cfgSingle.NumPaths) + 1j * randn(1, cfgSingle.NumPaths)) / sqrt(2 * cfgSingle.NumPaths);

            physicalChannelMatrix = LtvChannel(cfgSingle.TotalSubcarriers, pathDelays, dopplers, gains);

            % --- Single EP 链路 ---
            [txSigS, txDataS] = txSingle.transmit(pilotPowerSingle);
            noiseS = sqrt(noisePower / 2) * (randn(size(txSigS)) + 1j * randn(size(txSigS)));
            rxSigS = physicalChannelMatrix * txSigS + noiseS;
            rxDataS = rxSingle.receive(rxSigS, noisePower, physicalChannelMatrix, pilotPowerSingle);

            [eS, ~] = biterr(txDataS, rxDataS);
            errS = errS + eS;
            bitS = bitS + cfgSingle.NumActiveCarriers * cfgSingle.BitsPerSymbol;

            % --- CAZAC 链路 ---
            [txSigC, txDataC] = txCazac.transmit(pilotPowerCazac);
            noiseC = sqrt(noisePower / 2) * (randn(size(txSigC)) + 1j * randn(size(txSigC)));
            rxSigC = physicalChannelMatrix * txSigC + noiseC;
            rxDataC = rxCazac.receive(rxSigC, noisePower, physicalChannelMatrix, pilotPowerCazac);

            [eC, ~] = biterr(txDataC, rxDataC);
            errC = errC + eC;
            bitC = bitC + cfgCazac.NumActiveCarriers * cfgCazac.BitsPerSymbol;

            frameCount = frameCount + 1;

            if mod(frameCount, 500) == 0
                fprintf('  ZP=%2d, SNR=%2ddB | 帧 %5d | EP: %5d err | CAZAC: %5d err\r', ...
                    currentZp, snrDb, frameCount, errS, errC);
            end

        end

        berSingle(zpIdx, snrIdx) = errS / bitS;
        berCazac(zpIdx, snrIdx) = errC / bitC;

        fprintf('ZP=%2d, SNR=%2ddB | 帧数=%5d | EP BER=%.4e | CAZAC BER=%.4e\n', ...
            currentZp, snrDb, frameCount, berSingle(zpIdx, snrIdx), berCazac(zpIdx, snrIdx));
    end

    fprintf('\n');
end

%% ====================================================================
%  第四部分: 绘图
%  ====================================================================

colors = lines(numSnr);
markersSingle = {'o', 's', 'd'};
markersCazac = {'^', 'v', '>'};

% --- 图 1: BER vs ZP 长度 ---
figure('Position', [50 50 800 550]);
hold on; grid on;

for snrIdx = 1:numSnr
    semilogy(zpValues, berSingle(:, snrIdx), '-', ...
        'Color', colors(snrIdx, :), 'LineWidth', 1.8, ...
        'Marker', markersSingle{snrIdx}, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors(snrIdx, :));
    semilogy(zpValues, berCazac(:, snrIdx), '--', ...
        'Color', colors(snrIdx, :), 'LineWidth', 1.8, ...
        'Marker', markersCazac{snrIdx}, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'none');
end

xlabel('Zero-Padding 长度 (ZP)');
ylabel('BER');
title('缩减 ZP 对 BER 的影响: Single EP vs CAZAC + SIC');

legendEntries = cell(2 * numSnr, 1);

for snrIdx = 1:numSnr
    legendEntries{2 * snrIdx - 1} = sprintf('Single EP, SNR=%ddB', snrValues(snrIdx));
    legendEntries{2 * snrIdx} = sprintf('CAZAC+SIC, SNR=%ddB', snrValues(snrIdx));
end

legend(legendEntries, 'Location', 'northwest');
set(gca, 'YScale', 'log', 'XDir', 'reverse');
ylim([1e-5, 1]);
xline(77, ':', '理论完整 ZP', 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left');

% --- 图 2: BER vs 频谱效率 ---
figure('Position', [100 100 800 550]);
hold on; grid on;

for snrIdx = 1:numSnr
    seSingle = activeCarriersSingle ./ 256 * 2;
    seCazac = activeCarriersCazac ./ 256 * 2;

    semilogy(seSingle, berSingle(:, snrIdx), '-', ...
        'Color', colors(snrIdx, :), 'LineWidth', 1.8, ...
        'Marker', markersSingle{snrIdx}, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors(snrIdx, :));
    semilogy(seCazac, berCazac(:, snrIdx), '--', ...
        'Color', colors(snrIdx, :), 'LineWidth', 1.8, ...
        'Marker', markersCazac{snrIdx}, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'none');
end

xlabel('频谱效率 (bits/s/Hz, 归一化)');
ylabel('BER');
title('频谱效率 vs BER: Single EP vs CAZAC + SIC');
legend(legendEntries, 'Location', 'northeast');
set(gca, 'YScale', 'log');
ylim([1e-5, 1]);

%% ====================================================================
%  第五部分: 汇总表格
%  ====================================================================

fprintf('\n========================================================\n');
fprintf('  汇总: BER vs ZP 长度\n');
fprintf('========================================================\n');

for snrIdx = 1:numSnr
    fprintf('\n--- SNR = %d dB ---\n', snrValues(snrIdx));
    fprintf('%4s | %8s %8s | %10s %10s | %12s %12s\n', ...
        'ZP', 'EP开销', 'CA开销', 'EP有效载波', 'CA有效载波', 'EP BER', 'CA BER');
    fprintf(repmat('-', 1, 78)); fprintf('\n');

    for zpIdx = 1:numZp
        fprintf('%4d | %7.1f%% %7.1f%% | %10d %10d | %12.4e %12.4e\n', ...
            zpValues(zpIdx), ...
            overheadSingle(zpIdx) * 100, overheadCazac(zpIdx) * 100, ...
            activeCarriersSingle(zpIdx), activeCarriersCazac(zpIdx), ...
            berSingle(zpIdx, snrIdx), berCazac(zpIdx, snrIdx));
    end

end
