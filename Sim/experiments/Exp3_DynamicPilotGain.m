function resultStruct = Exp3_DynamicPilotGain( ...
        numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
        maxDopplerIdx, dopplerGuard, numPaths, pilotSnrDb, snrVecDb, figuresDir)
% EXP3_DYNAMICPILOTGAIN 执行固定导频与动态导频的全接收机对比实验。
%
%   描述:
%   在分数 Doppler 场景下保持 CLIP 接收机链路不变，仅改变导频功率策略，
%   对比固定导频与动态导频在 BER、NMSE 与 Pd 上的差异。
%
%   语法:
%   resultStruct = Exp3_DynamicPilotGain( ...
%       numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
%       maxDopplerIdx, dopplerGuard, numPaths, pilotSnrDb, snrVecDb, figuresDir)
%
%   输出:
%   resultStruct - (struct) 固定/动态导频的 BER、NMSE、Pd 与置信区间结果。
%
%   版本历史:
%   2026-04-22 - Aiden - 新增动态导频架构级增益实验。

    [fixedConfig, dynamicConfig] = createClipPilotConfigPair( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, ...
        dopplerGuard, numPaths, pilotSnrDb);
    fixedSystem = GiFreeSystem(fixedConfig);
    dynamicSystem = GiFreeSystem(dynamicConfig);

    trialCountVec = buildPilotGainTrialCountVec(snrVecDb);
    resultStruct = initializePilotGainResult(snrVecDb);

    fprintf('Exp3_DynamicPilotGain (P=%d, Fractional): 固定导频 vs 动态导频\n', numPaths);
    fprintf(['SNR   BER_Fixed   BER_Dynamic   NMSE_Fixed(dB)   NMSE_Dynamic(dB)' ...
        '   Pd_Fixed   Pd_Dynamic\n']);
    for snrIdx = 1:numel(snrVecDb)
        snrStats = runSinglePilotGainPoint( ...
            fixedSystem, dynamicSystem, bitsPerSymbol, snrVecDb(snrIdx), trialCountVec(snrIdx));
        resultStruct = writePilotGainPoint(resultStruct, snrIdx, snrStats);
        fprintf(' %2d   %.2e    %.2e       %+.2f            %+.2f          %.2f       %.2f\n', ...
            snrVecDb(snrIdx), snrStats.avgBer(1), snrStats.avgBer(2), ...
            10 * log10(snrStats.avgNmse(1) + 1e-20), ...
            10 * log10(snrStats.avgNmse(2) + 1e-20), ...
            snrStats.avgPd(1), snrStats.avgPd(2));
    end

    plotDynamicPilotGain(resultStruct, figuresDir);
end

function trialCountVec = buildPilotGainTrialCountVec(snrVecDb)
% BUILDPILOTGAINTRIALCOUNTVEC 返回动态导频实验的 SNR 自适应 trial 数量。

    trialCountVec = 300 * ones(numel(snrVecDb), 1);
    trialCountVec(snrVecDb <= 5) = 500;
    trialCountVec(snrVecDb >= 15) = 800;
end

function resultStruct = initializePilotGainResult(snrVecDb)
% INITIALIZEPILOTGAINRESULT 初始化动态导频实验结果结构体。

    numStrategies = 2;
    numPoints = numel(snrVecDb);
    resultStruct.snrVecDb = snrVecDb(:);
    resultStruct.strategyLabels = {'Fixed Pilot', 'Dynamic Pilot'};
    resultStruct.avgBer = zeros(numPoints, numStrategies);
    resultStruct.avgNmse = zeros(numPoints, numStrategies);
    resultStruct.avgPd = zeros(numPoints, numStrategies);
    resultStruct.berCiLower = zeros(numPoints, numStrategies);
    resultStruct.berCiUpper = zeros(numPoints, numStrategies);
end

function snrStats = runSinglePilotGainPoint( ...
        fixedSystem, dynamicSystem, bitsPerSymbol, dataSnrDb, trialCount)
% RUNSINGLEPILOTGAINPOINT 统计一个 SNR 点下两种导频策略的整机表现。

    noisePowerLin = 1;
    dataSnrLin = 10 ^ (dataSnrDb / 10);
    errorCount = zeros(1, 2);
    nmseAccum = zeros(1, 2);
    pdAccum = zeros(1, 2);
    totalBitCount = 0;

    for trialIdx = 1:trialCount
        trialRngState = rng;
        [truePathParams, trueEffectiveChannel] = dynamicSystem.ChannelSampler.sampleChannel();
        noiseVec = sqrt(noisePowerLin / 2) * ...
            (randn(dynamicSystem.Config.NumSubcarriers, 1) + ...
            1j * randn(dynamicSystem.Config.NumSubcarriers, 1));

        [txFrameFixed, txDataFixed] = transmitWithFixedSeed( ...
            fixedSystem, dataSnrLin, trialRngState);
        [txFrameDynamic, txDataDynamic] = transmitWithFixedSeed( ...
            dynamicSystem, dataSnrLin, trialRngState);
        validateSharedPayload(txDataFixed, txDataDynamic);

        rxSignalFixed = trueEffectiveChannel * txFrameFixed + noiseVec;
        rxSignalDynamic = trueEffectiveChannel * txFrameDynamic + noiseVec;

        [detectedFixed, nmseFixed, ~, diagFixed] = fixedSystem.Receiver.receive( ...
            rxSignalFixed, trueEffectiveChannel, dataSnrLin, noisePowerLin);
        [detectedDynamic, nmseDynamic, ~, diagDynamic] = dynamicSystem.Receiver.receive( ...
            rxSignalDynamic, trueEffectiveChannel, dataSnrLin, noisePowerLin);
        pdFixed = computePathDetectionProbability(truePathParams, diagFixed.finalEstimatedPaths);
        pdDynamic = computePathDetectionProbability(truePathParams, diagDynamic.finalEstimatedPaths);

        [errorCount(1), totalBitCount] = accumulatePilotGainErrors( ...
            txDataFixed, detectedFixed, bitsPerSymbol, errorCount(1), totalBitCount);
        [errorCount(2), ~] = accumulatePilotGainErrors( ...
            txDataDynamic, detectedDynamic, bitsPerSymbol, errorCount(2), 0);
        nmseAccum = nmseAccum + [nmseFixed, nmseDynamic];
        pdAccum = pdAccum + [pdFixed, pdDynamic];
    end

    snrStats.avgBer = [berFloor(errorCount(1), totalBitCount), ...
        berFloor(errorCount(2), totalBitCount)];
    snrStats.avgNmse = nmseAccum / trialCount;
    snrStats.avgPd = pdAccum / trialCount;
    for strategyIdx = 1:2
        [snrStats.berCiLower(strategyIdx), snrStats.berCiUpper(strategyIdx)] = ...
            computeConfidenceInterval(errorCount(strategyIdx), totalBitCount);
    end
end

function validateSharedPayload(txDataFixed, txDataDynamic)
% VALIDATESHAREDPAYLOAD 校验两种导频策略使用同一份数据符号索引。

    if ~isequal(txDataFixed, txDataDynamic)
        error('Exp3_DynamicPilotGain:PayloadMismatch', ...
            '固定导频与动态导频未共享同一组发送数据。');
    end
end

function [errorCount, totalBitCount] = accumulatePilotGainErrors( ...
        txDataIndices, detectedIndices, bitsPerSymbol, errorCount, totalBitCount)
% ACCUMULATEPILOTGAINERRORS 累积动态导频实验的误比特统计。

    [bitErrors, ~] = biterr(txDataIndices, detectedIndices, bitsPerSymbol);
    errorCount = errorCount + bitErrors;
    if totalBitCount > 0
        totalBitCount = totalBitCount + numel(txDataIndices) * bitsPerSymbol;
    else
        totalBitCount = numel(txDataIndices) * bitsPerSymbol;
    end
end

function resultStruct = writePilotGainPoint(resultStruct, snrIdx, snrStats)
% WRITEPILOTGAINPOINT 写回一个 SNR 点的动态导频实验结果。

    resultStruct.avgBer(snrIdx, :) = snrStats.avgBer;
    resultStruct.avgNmse(snrIdx, :) = snrStats.avgNmse;
    resultStruct.avgPd(snrIdx, :) = snrStats.avgPd;
    resultStruct.berCiLower(snrIdx, :) = snrStats.berCiLower;
    resultStruct.berCiUpper(snrIdx, :) = snrStats.berCiUpper;
end

function pdValue = computePathDetectionProbability(truePathParams, estimatedPaths)
% COMPUTEPATHDETECTIONPROBABILITY 计算估计路径相对真实路径的检测概率。

    if isempty(truePathParams)
        pdValue = 0;
        return;
    end

    if isempty(estimatedPaths)
        pdValue = 0;
        return;
    end

    matchedFlags = false(size(truePathParams, 1), 1);
    for estimatedIdx = 1:size(estimatedPaths, 1)
        for trueIdx = 1:size(truePathParams, 1)
            isSameDelay = estimatedPaths(estimatedIdx, 1) == truePathParams(trueIdx, 1);
            isSameDoppler = round(estimatedPaths(estimatedIdx, 2)) == round(truePathParams(trueIdx, 2));
            if ~matchedFlags(trueIdx) && isSameDelay && isSameDoppler
                matchedFlags(trueIdx) = true;
                break;
            end
        end
    end

    pdValue = sum(matchedFlags) / size(truePathParams, 1);
end

function plotDynamicPilotGain(resultStruct, figuresDir)
% PLOTDYNAMICPILOTGAIN 绘制动态导频增益实验图。

    ps = PlotStyle();
    colors = {ps.Blue, ps.Red};
    lineSpecs = {'--s', '-o'};
    figureHandle = figure('Position', [80 80 980 360], 'Color', 'w');

    axBer = subplot(1, 2, 1);
    hold(axBer, 'on');
    grid(axBer, 'on');
    box(axBer, 'on');
    for strategyIdx = 1:2
        ciLower = max(resultStruct.berCiLower(:, strategyIdx), 1e-6);
        ciUpper = max(resultStruct.berCiUpper(:, strategyIdx), 1e-6);
        patch(axBer, ...
            [resultStruct.snrVecDb.', fliplr(resultStruct.snrVecDb.')], ...
            [ciLower.', fliplr(ciUpper.')], colors{strategyIdx}, ...
            'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        semilogy(axBer, resultStruct.snrVecDb, max(resultStruct.avgBer(:, strategyIdx), 1e-6), ...
            lineSpecs{strategyIdx}, 'Color', colors{strategyIdx}, ...
            'LineWidth', 1.6, 'MarkerSize', 6, ...
            'DisplayName', resultStruct.strategyLabels{strategyIdx});
    end
    applyIeeeStyle(axBer, 9);
    set(axBer, 'YScale', 'log');
    ylim(axBer, [1e-6, 1]);
    xlabel(axBer, 'Data SNR (dB)');
    ylabel(axBer, 'BER');
    title(axBer, '(a) BER');
    legend(axBer, 'Location', 'southwest', 'FontSize', 8);

    axNmse = subplot(1, 2, 2);
    hold(axNmse, 'on');
    grid(axNmse, 'on');
    box(axNmse, 'on');
    for strategyIdx = 1:2
        plot(axNmse, resultStruct.snrVecDb, ...
            10 * log10(resultStruct.avgNmse(:, strategyIdx) + 1e-20), ...
            lineSpecs{strategyIdx}, 'Color', colors{strategyIdx}, ...
            'LineWidth', 1.6, 'MarkerSize', 6, ...
            'DisplayName', resultStruct.strategyLabels{strategyIdx});
    end
    applyIeeeStyle(axNmse, 9);
    xlabel(axNmse, 'Data SNR (dB)');
    ylabel(axNmse, 'NMSE (dB)');
    title(axNmse, '(b) Channel Estimation NMSE');
    legend(axNmse, 'Location', 'northeast', 'FontSize', 8);

    sgtitle('Fig.3: 动态导频的架构级增益', ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    savefig(figureHandle, fullfile(figuresDir, 'Fig3_DynamicPilotGain.fig'));
end
