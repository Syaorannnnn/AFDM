function resultStruct = Exp2_CfarThreshold( ...
        systemObj, modulationOrder, bitsPerSymbol, figuresDir, pathCountVecOverride)
% EXP2_CFARTHRESHOLD 执行阶段一检测器归因实验。
%
%   描述:
%   在固定导频前提下扫描 SNR × P × 场景 三维网格，对比四种检测方案：
%   Zhou 迭代+固定门限、Zhou 迭代+CFAR、OMP+固定门限、OMP+CFAR。
%   分数 Doppler 为主图，整数 Doppler 保留为不退化参考。
%
%   语法:
%   resultStruct = Exp2_CfarThreshold(systemObj, modulationOrder, bitsPerSymbol, figuresDir)
%
%   输入:
%   systemObj       - (GiFreeSystem) GI-Free 系统对象（整数 Doppler 配置）
%   modulationOrder - (double) 调制阶数
%   bitsPerSymbol   - (double) 每符号比特数
%   figuresDir      - (char/string) 图像输出目录
%   pathCountVecOverride - (double) 可选，定向运行的路径数向量
%
%   输出:
%   resultStruct - (struct) 包含两个场景（integer/fractional）、两个路径数（P=3/5）、
%                   四种方法的 BER/Pd/NMSE 统计及置信区间
%   resultStruct.snrVecDb      - (1x5 double) SNR 向量 [0, 5, 10, 15, 20]
%   resultStruct.pathCountVec  - (1x2 double) 路径数向量 [3, 5]
%   resultStruct.scenarioLabels - (1x2 cell) 场景标签 {'Integer', 'Fractional'}
%   resultStruct.methodLabels  - (1x4 cell) 方法标签
%   resultStruct.avgBer        - (5x2x2x4 double) BER 统计
%   resultStruct.avgPd         - (5x2x2x4 double) 检测概率统计
%   resultStruct.avgNmse       - (5x2x2x4 double) NMSE 统计
%   resultStruct.berCiLower    - (5x2x2x4 double) BER 置信区间下界
%   resultStruct.berCiUpper    - (5x2x2x4 double) BER 置信区间上界
%
%   版本历史:
%   2026-04-19 - Aiden - 新增 Part 2 实验函数。
%   2026-04-21 - Aiden - 重构为论文 Pt2 数据结构，增加分数 Doppler 场景。
%   2026-04-22 - Aiden - 调整为固定导频归因实验，并输出四方法完整曲线。

%% 参数设置
configObj = systemObj.Config;
snrVecDb = 0:5:20;
pathCountVec = resolveExp2PathCountVec(pathCountVecOverride);
scenarioLabels = {'Integer', 'Fractional'};
methodLabels = {'Zhou iterative + Fixed', 'Zhou iterative + CFAR', 'OMP + Fixed', 'OMP + CFAR'};
numMethods = 4;
numScenarios = 2;
numSnrPoints = numel(snrVecDb);
numPathCounts = numel(pathCountVec);

% 自适应 trial 数量
trialCountVec = buildExp2TrialCountVec(snrVecDb.');

%% 结果结构初始化
resultStruct.snrVecDb = snrVecDb;
resultStruct.pathCountVec = pathCountVec;
resultStruct.scenarioLabels = scenarioLabels;
resultStruct.methodLabels = methodLabels;
resultStruct.avgBer = zeros(numSnrPoints, numPathCounts, numScenarios, numMethods);
resultStruct.avgPd = zeros(numSnrPoints, numPathCounts, numScenarios, numMethods);
resultStruct.avgNmse = zeros(numSnrPoints, numPathCounts, numScenarios, numMethods);
resultStruct.berCiLower = zeros(numSnrPoints, numPathCounts, numScenarios, numMethods);
resultStruct.berCiUpper = zeros(numSnrPoints, numPathCounts, numScenarios, numMethods);

%% 主扫描循环
fprintf('Exp2_CfarThreshold: 固定导频阶段一检测器归因实验\n');
for pathCountIdx = 1:numPathCounts
    numPaths = pathCountVec(pathCountIdx);
    configObj.NumPaths = numPaths;
    configObj.NumPathsUpper = resolveExp2NumPathsUpper(numPaths);
    for scenarioIdx = 1:numScenarios
        isFractional = (scenarioIdx == 2);

        % 动态切换配置
        configObj.UseFractionalDoppler = isFractional;
        if isFractional
            configObj.DirichletRadius = 1;
        else
            configObj.DirichletRadius = 0;
        end

        scenarioName = scenarioLabels{scenarioIdx};
        fprintf('  P=%d, Scenario=%s\n', numPaths, scenarioName);

        for snrIdx = 1:numSnrPoints
            dataSnrDb = snrVecDb(snrIdx);
            trialCount = trialCountVec(snrIdx);

            pointStats = runSingleCfarPoint( ...
                systemObj, modulationOrder, bitsPerSymbol, dataSnrDb, trialCount, isFractional, numPaths);

            resultStruct.avgBer(snrIdx, pathCountIdx, scenarioIdx, :) = pointStats.avgBer;
            resultStruct.avgPd(snrIdx, pathCountIdx, scenarioIdx, :) = pointStats.avgPd;
            resultStruct.avgNmse(snrIdx, pathCountIdx, scenarioIdx, :) = pointStats.avgNmse;
            resultStruct.berCiLower(snrIdx, pathCountIdx, scenarioIdx, :) = pointStats.berCiLower;
            resultStruct.berCiUpper(snrIdx, pathCountIdx, scenarioIdx, :) = pointStats.berCiUpper;

            fprintf('    SNR=%2d dB: BER=[%.2e %.2e %.2e %.2e], Pd=[%.2f %.2f %.2f %.2f]\n', ...
                dataSnrDb, pointStats.avgBer, pointStats.avgPd);
        end
    end
end

%% 绘图
plotCfarBenchmark(resultStruct, figuresDir, configObj);

end

%% ==================== 单点统计函数 ====================

function pointStats = runSingleCfarPoint( ...
        systemObj, modulationOrder, bitsPerSymbol, dataSnrDb, trialCount, isFractional, numPaths)
% RUNSINGLECFARPOINT 统计一个 SNR×P×场景 网格点的四种方法表现。

configObj = systemObj.Config;
numSubcarriers = configObj.NumSubcarriers;
dataPos1 = configObj.DataPos1;

dataSnrLin = 10 ^ (dataSnrDb / 10);
[pilotFrame, pilotAmplitude] = buildGiFreePilotFrameForSnr(configObj, dataSnrLin);

% 固定门限计算：仅基于噪声功率（常数 3σ），不随 data SNR 变化
% 动态导频模式下 pilotAmplitude 已随 SNR 调整，门限应同步调整但不额外乘 dataSnrLin
paperAmpThreshold = 3;  % 3 倍噪声标准差（噪声功率归一化为 1）
paperMetricThreshold = (pilotAmplitude^2) * paperAmpThreshold^2;

% OMP 正则化参数
ompRegParam = (configObj.NumDataSymbols / numSubcarriers) * dataSnrLin / pilotAmplitude^2;

% LMMSE 正则化参数
lmmseRegParam = 1 / dataSnrLin;

% 统计量累积
accPd = zeros(1, 4);
accNmse = zeros(1, 4);
bitErrorCount = zeros(1, 4);
totalBitCount = 0;

%% Trial 循环
for trialIdx = 1:trialCount
    % 生成随机信道（共享给四种方法）
    [trueDelays, trueDopplers, trueGains] = generateRandomChannel( ...
        numPaths, configObj.MaxDelaySamples, configObj.MaxDopplerIdx, isFractional);
    truePathParams = [trueDelays, trueDopplers, trueGains];

    % 构建有效信道
    trueEffectiveChannel = systemObj.ChannelBuilder.buildEffectiveChannel(truePathParams);
    trueNorm = max(norm(full(trueEffectiveChannel), 'fro')^2, 1e-20);

    % 发射与接收
    [txFrame, txDataIndices] = systemObj.Transmitter.transmit(dataSnrLin);
    noiseVec = sqrt(0.5) * (randn(numSubcarriers, 1) + 1j * randn(numSubcarriers, 1));
    rxSignal = trueEffectiveChannel * txFrame + noiseVec;

    totalBitCount = totalBitCount + numel(txDataIndices) * bitsPerSymbol;

    % CFAR 门限计算（基于当前接收信号）
    cfarMetricThreshold = paperCfarThreshold( ...
        rxSignal, numSubcarriers, pilotAmplitude, ...
        configObj.MaxDelaySamples, configObj.MaxDopplerIdx, ...
        configObj.CfarFalseAlarmProb);

    % 执行四种检测方法
    [detectedCell, channelCell] = runFourThresholdMethods( ...
        systemObj, rxSignal, paperMetricThreshold, cfarMetricThreshold, ompRegParam, configObj.NumPathsUpper);

    % 累积统计量
    numTruePaths = size(truePathParams, 1);
    for methodIdx = 1:4
        [numCorrect, ~] = matchPathsLocal(truePathParams, detectedCell{methodIdx});
        accPd(methodIdx) = accPd(methodIdx) + numCorrect / numTruePaths;
        accNmse(methodIdx) = accNmse(methodIdx) + ...
            norm(full(channelCell{methodIdx}) - trueEffectiveChannel, 'fro')^2 / trueNorm;

        bitErrorCount(methodIdx) = bitErrorCount(methodIdx) + quickBitErrorsLocal( ...
            rxSignal, channelCell{methodIdx}, pilotFrame, lmmseRegParam, ...
            numSubcarriers, dataPos1, dataSnrLin, modulationOrder, txDataIndices);
    end
end

%% 输出统计结果
pointStats.avgPd = accPd / trialCount;
pointStats.avgNmse = accNmse / trialCount;
pointStats.avgBer = bitErrorCount / totalBitCount;
pointStats.berCiLower = zeros(1, 4);
pointStats.berCiUpper = zeros(1, 4);

for methodIdx = 1:4
    [pointStats.berCiLower(methodIdx), pointStats.berCiUpper(methodIdx)] = ...
        computeConfidenceInterval(bitErrorCount(methodIdx), totalBitCount);
end

end

%% ==================== 四种方法执行 ====================

function [detectedCell, channelCell] = runFourThresholdMethods( ...
        systemObj, rxSignal, paperMetricThreshold, cfarMetricThreshold, ompRegParam, maxPaths)
% RUNFOURTHRESHOLDMETHODS 执行四种门限方案。

configObj = systemObj.Config;
channelBuilder = systemObj.ChannelBuilder;
estimator = systemObj.Estimator;

detectedCell = cell(4, 1);
channelCell = cell(4, 1);

% 方法 1: Zhou 迭代 + 固定门限
[detectedCell{1}, channelCell{1}] = zhouIterativeDetect( ...
    rxSignal, configObj, channelBuilder, paperMetricThreshold, maxPaths);

% 方法 2: Zhou 迭代 + CFAR 门限
[detectedCell{2}, channelCell{2}] = zhouIterativeDetect( ...
    rxSignal, configObj, channelBuilder, cfarMetricThreshold, maxPaths);

% 方法 3: OMP + 固定门限
weights = ones(configObj.NumSubcarriers, 1);
[detectedCell{3}, channelCell{3}] = estimator.estimateByOmpWithFixedThreshold( ...
    rxSignal, ompRegParam, weights, paperMetricThreshold);

% 方法 4: OMP + CFAR（内部自适应门限）
[detectedCell{4}, channelCell{4}] = estimator.estimateByOmp(rxSignal, ompRegParam, weights);

end

%% ==================== 辅助函数 ====================

function [numCorrect, numFalse] = matchPathsLocal(truePaths, estimatedPaths)
% MATCHPATHSLOCAL 按 delay 与 round(doppler) 做路径匹配。

if isempty(truePaths)
    numCorrect = 0;
    numFalse = size(estimatedPaths, 1);
    return;
end

if isempty(estimatedPaths)
    numCorrect = 0;
    numFalse = 0;
    return;
end

matchedFlags = false(size(truePaths, 1), 1);
falseAlarmFlags = true(size(estimatedPaths, 1), 1);

for estimatedIdx = 1:size(estimatedPaths, 1)
    for trueIdx = 1:size(truePaths, 1)
        isSameDelay = estimatedPaths(estimatedIdx, 1) == truePaths(trueIdx, 1);
        isSameDoppler = round(estimatedPaths(estimatedIdx, 2)) == round(truePaths(trueIdx, 2));
        if ~matchedFlags(trueIdx) && isSameDelay && isSameDoppler
            matchedFlags(trueIdx) = true;
            falseAlarmFlags(estimatedIdx) = false;
            break;
        end
    end
end

numCorrect = sum(matchedFlags);
numFalse = sum(falseAlarmFlags);

end

function numBitErrors = quickBitErrorsLocal( ...
        rxSignal, estimatedChannel, pilotFrame, regParam, numSubcarriers, ...
        dataPos1, dataSnrLin, modulationOrder, txDataIndices)
% QUICKBITERRORSLOCAL 计算一次 LMMSE 检测的误比特数。

cleanSignal = rxSignal - estimatedChannel * pilotFrame;
estimatedSignal = GiFreeReceiver.lmmseDetect( ...
    cleanSignal, estimatedChannel, regParam, numSubcarriers, dataPos1);
detectedIndices = qamdemod( ...
    estimatedSignal(dataPos1) / sqrt(dataSnrLin), modulationOrder, 'UnitAveragePower', true);
[numBitErrors, ~] = biterr(txDataIndices, detectedIndices, log2(modulationOrder));

end

%% ==================== 绘图函数 ====================

function plotCfarBenchmark(resultStruct, figuresDir, ~)
% PLOTCFARBENCHMARK 绘制阶段一检测器归因实验图。

styleList = buildPhase1StyleList();
scenarioFileSuffix = {'Integer', 'Fractional'};
minBerValue = 1e-6;

for scenarioIdx = 1:2
    figureHandle = figure('Position', [80 80 980 360], 'Color', 'w');
    for pathCountIdx = 1:numel(resultStruct.pathCountVec)
        axHandle = subplot(1, numel(resultStruct.pathCountVec), pathCountIdx);
        plotPhase1Panel(axHandle, resultStruct, scenarioIdx, pathCountIdx, styleList, minBerValue);
    end
    sgtitle(sprintf('Fig.2: 阶段一检测器归因实验（%s Doppler）', ...
        resultStruct.scenarioLabels{scenarioIdx}), ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    savefig(figureHandle, fullfile(figuresDir, ...
        sprintf('Fig2_Phase1Ablation_%s.fig', scenarioFileSuffix{scenarioIdx})));
end

fprintf(['  Part 2 figures saved: Fig2_Phase1Ablation_Integer.fig, ', ...
    'Fig2_Phase1Ablation_Fractional.fig\n']);

end

function styleList = buildPhase1StyleList()
% BUILDPHASE1STYLELIST 返回阶段一归因实验的曲线样式。

    ps = PlotStyle();
    styleList = { ...
        struct('Color', ps.Blue,   'LineSpec', '--s'), ...
        struct('Color', ps.Purple, 'LineSpec', '-d'), ...
        struct('Color', ps.Green,  'LineSpec', '-.^'), ...
        struct('Color', ps.Red,    'LineSpec', '-o')};
end

function plotPhase1Panel(axHandle, resultStruct, scenarioIdx, pathCountIdx, styleList, minBerValue)
% PLOTPHASE1PANEL 绘制单个场景与路径数的四方法 BER 曲线。

    hold(axHandle, 'on');
    grid(axHandle, 'on');
    box(axHandle, 'on');

    snrVecDb = resultStruct.snrVecDb;
    for methodIdx = 1:numel(resultStruct.methodLabels)
        berValues = max(squeeze(resultStruct.avgBer(:, pathCountIdx, scenarioIdx, methodIdx)), minBerValue);
        ciLower = max(squeeze(resultStruct.berCiLower(:, pathCountIdx, scenarioIdx, methodIdx)), minBerValue);
        ciUpper = max(squeeze(resultStruct.berCiUpper(:, pathCountIdx, scenarioIdx, methodIdx)), minBerValue);
        styleObj = styleList{methodIdx};

        patch(axHandle, ...
            [snrVecDb, fliplr(snrVecDb)], ...
            [ciLower.', fliplr(ciUpper.')], styleObj.Color, ...
            'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        semilogy(axHandle, snrVecDb, berValues, styleObj.LineSpec, ...
            'Color', styleObj.Color, 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'DisplayName', resultStruct.methodLabels{methodIdx});
    end

    applyIeeeStyle(axHandle, 9);
    set(axHandle, 'YScale', 'log');
    ylim(axHandle, [minBerValue, 1]);
    xlabel(axHandle, 'Data SNR (dB)');
    ylabel(axHandle, 'BER');
    title(axHandle, sprintf('P=%d', resultStruct.pathCountVec(pathCountIdx)));
    legend(axHandle, 'Location', 'southwest', 'FontSize', 8);
end

