% MainGiFreeSim.m — GI-Free AFDM 最终系统仿真 (handle 类架构)
%
% 对比两条核心曲线:
%   Perfect CSI    —— 已知真实信道的理论性能上界
%   GiFreeSystem   —— 最终版 Turbo DD-CE 接收机
%
% ============================================================
%   与静态方法版本的核心区别:
%     1. 所有模块通过 handle 引用共享同一个 GiFreeConfig
%     2. GiFreeSystem 在构造时一站式组装所有子模块
%     3. 运行时只需 sys.runTrial(...), 无需传递 cfg
%     4. 修改 cfg 属性后, 所有模块自动感知新值
% ============================================================
%
% 图表配色: 深色莫兰迪 (Morandi) 色系, 纯白背景

clear; clc; close all; rng(42);
addpath(genpath('./'));
%% ======================== 1. 仿真参数 ========================

% ---- 构造共享配置 (handle 对象, 全系统唯一实例) ----
cfg = GiFreeConfig();
cfg.NumSubcarriers       = 2^9;
cfg.ModulationOrder      = 4;
cfg.MaxDelaySpread       = 2;
cfg.MaxDopplerIndex      = 2;
cfg.NumPaths             = 3;
cfg.PilotSnrDb           = 45;
cfg.MaxSicIterations     = 10;  % Phase1(3) + Phase2(6) + PostDecision(1)
cfg.NumPathsUpper        = cfg.NumPaths + 1;
cfg.DopplerGuard         = 4;
cfg.SpreadWidth          = 4;
cfg.UseFractionalDoppler = true;

% ---- 一站式构造完整系统 ----
%   GiFreeSystem 构造函数会自动创建:
%     Transmitter, ChannelBuilder, Estimator, Receiver
%   并将它们全部连接到同一个 cfg handle.
sys = GiFreeSystem(cfg);

% ---- 仿真控制参数 ----
snrDbList  = [0, 5, 10, 15, 20];
numTrials  = 200;
noisePower = 1;

numSnrPts  = length(snrDbList);
berResults = zeros(numSnrPts, 2);
mseResults = zeros(numSnrPts, 2);

fprintf('============================================================\n');
fprintf('  GI-Free AFDM — Final System Simulation (Handle Architecture)\n');
fprintf('  N = %d, %d-QAM, P = %d, PilotSNR = %ddB, SIC = %d iter\n', ...
    cfg.NumSubcarriers, cfg.ModulationOrder, cfg.NumPaths, ...
    cfg.PilotSnrDb, cfg.MaxSicIterations);
fprintf('============================================================\n');

%% ======================== 2. 蒙特卡洛主循环 ========================

for snrIdx = 1:numSnrPts
    currentSnrDb  = snrDbList(snrIdx);
    dataSnrLinear = 10 ^ (currentSnrDb / 10);
    bitsPerSymbol = log2(cfg.ModulationOrder);

    cumBitErrSys = 0;
    cumBitErrRef = 0;
    cumTotalBits = 0;
    cumMse       = 0;

    for trialIdx = 1:numTrials
        % 核心调用: 实例方法, 无需传递 cfg
        result = sys.runTrial(dataSnrLinear, noisePower);

        cumBitErrSys = cumBitErrSys + result.bitErrorsSys;
        cumBitErrRef = cumBitErrRef + result.bitErrorsRef;
        cumTotalBits = cumTotalBits + result.totalBits;
        cumMse       = cumMse + result.mseSystem;
    end

    % 统计置信下限: 0 个错误时用 1/(2*totalBits) 代替
    berFloor = 1 / (2 * max(cumTotalBits, 1));

    rawBerRef = cumBitErrRef / max(cumTotalBits, 1);
    rawBerSys = cumBitErrSys / max(cumTotalBits, 1);

    berResults(snrIdx, 1) = max(rawBerRef, berFloor);
    berResults(snrIdx, 2) = max(rawBerSys, berFloor);
    mseResults(snrIdx, 2) = cumMse / numTrials;

    % 标记触及仿真精度边界的点
    isFlooredRef(snrIdx) = (cumBitErrRef == 0);  %#ok<SAGROW>
    isFlooredSys(snrIdx) = (cumBitErrSys == 0);  %#ok<SAGROW>

    berGap = berResults(snrIdx, 2) / berResults(snrIdx, 1);
    if isFlooredRef(snrIdx)
        fprintf('SNR = %2d dB | BER: PerfCSI < %.2e* System = %.2e  (Gap N/A)  | NMSE = %.2e\n', ...
            currentSnrDb, berFloor, berResults(snrIdx, 2), mseResults(snrIdx, 2));
    else
        fprintf('SNR = %2d dB | BER: PerfCSI = %.2e  System = %.2e  (Gap = %.1f x) | NMSE = %.2e\n', ...
            currentSnrDb, berResults(snrIdx, 1), berResults(snrIdx, 2), ...
            berGap, mseResults(snrIdx, 2));
    end
end

%% ======================== 3. 结果可视化 ========================

plotFinalResults(snrDbList, berResults, mseResults, cfg, isFlooredRef, isFlooredSys);
savefig(gcf, 'GiFreeFinalResults.fig');

%% ======================== 辅助函数: 绘图 ========================

function plotFinalResults(snrDbList, berResults, mseResults, cfg, isFlooredRef, isFlooredSys)

    colPerfCsi   = [0.32, 0.30, 0.28];
    colSystem    = [0.33, 0.47, 0.38];
    colAccent    = [0.55, 0.38, 0.40];
    colGridLine  = [0.82, 0.80, 0.78];
    colFillArea  = [0.85, 0.90, 0.86];
    colBarFace   = [0.40, 0.55, 0.45];
    colBarEdge   = [0.28, 0.40, 0.32];

    fontLabel = 11;
    fontTitle = 12;
    fontLeg   = 9;

    numSc    = cfg.NumSubcarriers;
    modOrder = cfg.ModulationOrder;
    numP     = cfg.NumPaths;
    dopGuard = cfg.DopplerGuard;
    spreadKv = cfg.SpreadWidth;

    figure('Position', [80, 60, 1200, 900], 'Color', [1, 1, 1], ...
        'Name', 'GI-Free AFDM — Final System Performance');

    % ---- 子图 1: BER 核心对比 ----
    subplot(2, 2, 1);
    set(gca, 'Color', [1, 1, 1]);

    semilogy(snrDbList, berResults(:, 1), '-^', 'Color', colPerfCsi, ...
        'LineWidth', 2.2, 'MarkerSize', 9, 'MarkerFaceColor', colPerfCsi, ...
        'DisplayName', 'Perfect CSI'); hold on;
    semilogy(snrDbList, berResults(:, 2), '-s', 'Color', colSystem, ...
        'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
        'DisplayName', 'GiFreeSystem');

    % 标注仿真精度边界点: 空心标记 + 向下箭头
    for i = 1:length(snrDbList)
        if isFlooredRef(i)
            semilogy(snrDbList(i), berResults(i, 1), '^', 'Color', colPerfCsi, ...
                'MarkerSize', 11, 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, ...
                'HandleVisibility', 'off');
            text(snrDbList(i), berResults(i, 1) * 2.5, '\downarrow 0 errors', ...
                'HorizontalAlignment', 'center', 'FontSize', 7, ...
                'Color', colPerfCsi, 'FontWeight', 'bold');
        end
        if isFlooredSys(i)
            semilogy(snrDbList(i), berResults(i, 2), 's', 'Color', colSystem, ...
                'MarkerSize', 12, 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, ...
                'HandleVisibility', 'off');
            text(snrDbList(i), berResults(i, 2) * 2.5, '\downarrow 0 errors', ...
                'HorizontalAlignment', 'center', 'FontSize', 7, ...
                'Color', colSystem, 'FontWeight', 'bold');
        end
    end

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6, 'MinorGridAlpha', 0.3);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('BER', 'FontSize', fontLabel);
    title(sprintf('BER — N=%d, %d-QAM, P=%d', numSc, modOrder, numP), ...
        'FontSize', fontTitle);
    legend('Location', 'southwest', 'FontSize', fontLeg, 'Box', 'off');
    hold off;

    % ---- 子图 2: NMSE ----
    subplot(2, 2, 2);
    set(gca, 'Color', [1, 1, 1]);

    semilogy(snrDbList, mseResults(:, 2), '-s', 'Color', colSystem, ...
        'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
        'DisplayName', 'GiFreeSystem');

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('Normalized MSE', 'FontSize', fontLabel);
    title('Channel Estimation NMSE', 'FontSize', fontTitle);
    legend('Location', 'best', 'FontSize', fontLeg, 'Box', 'off');

    % ---- 子图 3: BER 差距分析 ----
    %   仅对两条曲线均有实测错误的 SNR 点计算有效 Gap;
    %   触及仿真边界的点标注 "N/A" 而非绘制失真的数值.
    subplot(2, 2, 3);
    set(gca, 'Color', [1, 1, 1]);

    berGap    = berResults(:, 2) ./ berResults(:, 1);
    validMask = ~isFlooredRef(:) & ~isFlooredSys(:);

    if any(validMask)
        validSnrs = snrDbList(validMask);
        validGap  = berGap(validMask);

        fill([validSnrs, fliplr(validSnrs)], ...
            [validGap', ones(1, sum(validMask))], ...
            colFillArea, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off'); hold on;
        plot(validSnrs, validGap, '-s', 'Color', colSystem, ...
            'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
            'DisplayName', 'BER_{System} / BER_{PerfCSI}');
    else
        hold on;
    end

    yline(1, ':', 'Color', colGridLine, 'LineWidth', 1.8, ...
        'DisplayName', 'Perfect CSI = 1');

    % 对触及边界的点做标注
    for i = 1:length(snrDbList)
        if isFlooredRef(i) || isFlooredSys(i)
            yPos = 1.1;
            text(snrDbList(i), yPos, 'N/A', ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', colAccent, 'FontWeight', 'bold');
        end
    end

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('BER / BER_{PerfCSI}', 'FontSize', fontLabel);
    title('BER Gap to Perfect CSI', 'FontSize', fontTitle);
    legend('Location', 'best', 'FontSize', fontLeg, 'Box', 'off');
    hold off;

    % ---- 子图 4: 选定 SNR 柱状图 ----
    subplot(2, 2, 4);
    set(gca, 'Color', [1, 1, 1]);

    selectedSnrs = [5, 10, 15, 20, 25];
    selIdx = zeros(size(selectedSnrs));
    for i = 1:length(selectedSnrs)
        [~, selIdx(i)] = min(abs(snrDbList - selectedSnrs(i)));
    end

    berForBar = zeros(length(selectedSnrs), 2);
    for i = 1:length(selectedSnrs)
        berForBar(i, :) = berResults(selIdx(i), :);
    end

    barHandle = bar(1:length(selectedSnrs), berForBar, 'grouped');
    set(gca, 'YScale', 'log');
    barHandle(1).FaceColor = colPerfCsi;
    barHandle(1).EdgeColor = colPerfCsi;
    barHandle(2).FaceColor = colBarFace;
    barHandle(2).EdgeColor = colBarEdge;

    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d dB', x), ...
        selectedSnrs, 'UniformOutput', false));

    for i = 1:length(selectedSnrs)
        si = selIdx(i);
        isAnyFloored = isFlooredRef(si) || isFlooredSys(si);

        if isAnyFloored
            % 触及仿真边界: 标注而非显示失真倍数
            barCenterX = i;
            barTopY    = max(berForBar(i, :)) * 2;
            text(barCenterX, barTopY, '\leq bound', ...
                'HorizontalAlignment', 'center', 'FontSize', 7, ...
                'FontWeight', 'bold', 'Color', colAccent);
        else
            gapVal = berForBar(i, 2) / berForBar(i, 1);
            if gapVal < 100 && berForBar(i, 2) > 0
                barCenterX = i + 0.15;
                barTopY    = berForBar(i, 2) * 1.5;
                text(barCenterX, barTopY, sprintf('%.1fx', gapVal), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, ...
                    'FontWeight', 'bold', 'Color', colAccent);
            end
        end
    end

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR', 'FontSize', fontLabel);
    ylabel('BER', 'FontSize', fontLabel);
    title('BER at Selected SNR Points', 'FontSize', fontTitle);
    legend('Perfect CSI', 'GiFreeSystem', ...
        'Location', 'northeast', 'FontSize', fontLeg, 'Box', 'off');

    sgtitle(sprintf( ...
        'GI-Free AFDM Final System — \\xi_\\nu=%d, k_\\nu=%d, Turbo DD-CE', ...
        dopGuard, spreadKv), 'FontSize', 14, 'FontWeight', 'bold', 'Color', colPerfCsi);
end