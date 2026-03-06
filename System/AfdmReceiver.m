classdef AfdmReceiver < handle
    % AfdmReceiver: AFDM 接收机，包含信道估计、均衡与干扰消除

    properties (Access = public)
        % --- 接收机配置 ---
        CsiMode (1, 1) string = "Estimated" % "Perfect" 或 "Estimated"
        EqualizerType (1, 1) string = "MMSE" % "MMSE", "MRC", 或 "MRC-DFE"

        % --- 迭代与 SIC 配置 ---
        EnableIterativeIc (1, 1) logical = true % 是否开启缩减 ZP 的迭代干扰消除
        NumIcIterations (1, 1) double = 3 % SIC 迭代次数
        NumMaxIterations (1, 1) double = 15 % MRC-DFE 均衡器最大迭代次数
    end

    properties (Access = private)
        Config % 全局配置对象引用
        currentPilotPower (1, 1) double = 1;
    end

    methods (Access = public)

        % --- 构造函数 ---
        function obj = AfdmReceiver(configObj)
            obj.Config = configObj;
        end

        % --- 接收主函数 ---
        function rxData = receive(obj, rxSignal, noisePower, physicalChannelMatrix, pilotPower)
            obj.currentPilotPower = pilotPower; % 记录当前功率
            % 移除前缀
            rxNoPrefix = rxSignal(obj.Config.PrefixLength + 1:end);

            % 变换到 DAFT / DFT 域
            if upper(obj.Config.WaveformType) == "AFDM"
                daftRx = AfdmTransforms.daft(rxNoPrefix, obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                daftRx = AfdmTransforms.dft(rxNoPrefix);
            end

            % 匹配滤波
            matchedDaftRx = daftRx .* conj(obj.Config.PulseShapingWindow);

            rxData = obj.processDetection(matchedDaftRx, noisePower, physicalChannelMatrix);
        end

    end

    methods (Access = private)

        function rxData = processDetection(obj, daftRx, noisePower, physicalChannelMatrix)

            rxDataDemod = zeros(length(obj.Config.ActiveIndices), 1);
            residualSignalForPilot = daftRx;
            qamConstellation = qammod(0:obj.Config.ModulationOrder - 1, obj.Config.ModulationOrder, 'UnitAveragePower', true);
            maxAmp = max(abs(qamConstellation));

            % 判断是否满足启用迭代 SIC 的条件
            isSicEnabled = (obj.CsiMode == "Estimated") && (upper(obj.Config.PilotType) == "CAZAC") && obj.EnableIterativeIc;

            % 动态设定迭代次数
            if isSicEnabled
                actualIters = obj.NumIcIterations;
            else
                actualIters = 1; % 理想 CSI 或标准模式只需跑一次
            end

            % 核心处理循环
            for iter = 1:actualIters

                % --- 1. 信道获取 ---
                if obj.CsiMode == "Perfect"
                    effectiveChannelMatrix = AfdmTransforms.generateEffectiveChannelMatrix(physicalChannelMatrix, obj.Config);
                else
                    [effectiveChannelMatrix, ~] = obj.runChannelEstimator(residualSignalForPilot);
                end

                % --- 消除IP2D ---
                if isSicEnabled
                    mockPilotFrame = zeros(obj.Config.NumDataSubcarriers, 1);
                    startIdx = obj.Config.PilotIndex;
                    endIdx = startIdx + obj.Config.PilotSequenceLength - 1;
                    mockPilotFrame(startIdx:endIdx) = obj.getScaledPilot();

                    ip2dInterference = effectiveChannelMatrix * mockPilotFrame;
                    residualSignalForData = daftRx - ip2dInterference;
                else
                    residualSignalForData = daftRx; % 非 SIC 模式直接送入均衡器
                end

                % --- 3. 数据均衡 ---
                eqSignal = obj.runEqualizer(residualSignalForData, effectiveChannelMatrix, noisePower);

                % --- 4. 判决与跳出 ---
                if iter == actualIters
                    % 最后一轮或唯一一轮：直接硬判决并输出
                    rxDataDemod = qamdemod(eqSignal, obj.Config.ModulationOrder, 'UnitAveragePower', true);
                    break;
                end

                % --- 迭代重构 ---
                if iter == 1
                    % 第一轮：软判决 + 限幅防发散
                    alpha = 0.6;
                    reModData = eqSignal;
                    exceedIdx = abs(reModData) > maxAmp;
                    reModData(exceedIdx) = maxAmp * exp(1j * angle(reModData(exceedIdx)));
                else
                    % 中间轮次：硬判决重构
                    alpha = 1.0;
                    tempDemod = qamdemod(eqSignal, obj.Config.ModulationOrder, 'UnitAveragePower', true);
                    reModData = qammod(tempDemod, obj.Config.ModulationOrder, 'UnitAveragePower', true);
                end

                % 消除ID2P，更新信道估计残差
                mockDataFrame = zeros(obj.Config.NumDataSubcarriers, 1);
                mockDataFrame(obj.Config.ActiveIndices) = reModData(:);
                id2pInterference = effectiveChannelMatrix * mockDataFrame;

                residualSignalForPilot = daftRx - alpha * id2pInterference;
            end

            rxData = rxDataDemod(:);
        end

        function eqSignal = runEqualizer(obj, rxSignal, effectiveChannelMatrix, noisePower)

            switch upper(obj.EqualizerType)
                case "MMSE"
                    eqSignal = obj.runMmse(rxSignal, effectiveChannelMatrix, noisePower);
                case {"MRC", "MRC-DFE"}
                    eqSignal = obj.runWeightedMrcDfe(rxSignal, effectiveChannelMatrix, noisePower);
                otherwise
                    error('AfdmReceiver:UnknownEqualizer', '未知的均衡器类型: %s', obj.EqualizerType);
            end

        end

        function estData = runMmse(obj, rxSignal, effectiveChannelMatrix, noisePower)
            activeChannelMatrix = effectiveChannelMatrix(:, obj.Config.ActiveIndices);
            gramMatrix = activeChannelMatrix' * activeChannelMatrix;
            estData = (gramMatrix + noisePower * eye(obj.Config.NumActiveCarriers)) \ (activeChannelMatrix' * rxSignal);
        end

        function estData = runWeightedMrcDfe(obj, rxSignal, effectiveChannelMatrix, noisePower)
            numSubcarriers = obj.Config.NumDataSubcarriers;
            epsilon = 1e-5;
            columnEnergies = full(sum(abs(effectiveChannelMatrix) .^ 2, 1)).';

            columnRowIndices = cell(numSubcarriers, 1);
            columnValues = cell(numSubcarriers, 1);

            for k = 1:numSubcarriers
                [rows, ~, values] = find(effectiveChannelMatrix(:, k));
                columnRowIndices{k} = rows;
                columnValues{k} = values;
            end

            currentEstimate = zeros(numSubcarriers, 1);
            previousEstimate = zeros(numSubcarriers, 1);
            residualSignal = rxSignal;

            activeEnergies = columnEnergies(obj.Config.ActiveIndices);
            [~, sortIdx] = sort(activeEnergies, 'descend');
            orderedIndices = obj.Config.ActiveIndices(sortIdx);

            for n = 1:obj.NumMaxIterations

                for k = orderedIndices
                    rows = columnRowIndices{k};
                    if isempty(rows), continue; end
                    channelValues = columnValues{k};

                    mrcOutput = sum(conj(channelValues) .* residualSignal(rows)) + columnEnergies(k) * previousEstimate(k);
                    newEstimate = mrcOutput / (columnEnergies(k) + noisePower);
                    estimateChange = newEstimate - previousEstimate(k);
                    currentEstimate(k) = newEstimate;

                    if abs(estimateChange) > 1e-6
                        residualSignal(rows) = residualSignal(rows) - channelValues * estimateChange;
                    end

                end

                if norm(currentEstimate - previousEstimate) < epsilon
                    break;
                else
                    previousEstimate = currentEstimate;
                end

            end

            estData = currentEstimate(obj.Config.ActiveIndices);
        end

        function [effectiveChannelMatrix, finalPaths] = runChannelEstimator(obj, rxSignalDaft)
            % 构造仅含导频的参考帧
            mockTxFrame = zeros(obj.Config.NumDataSubcarriers, 1);
            startIdx = obj.Config.PilotIndex;
            endIdx = startIdx + obj.Config.PilotSequenceLength - 1;
            mockTxFrame(startIdx:endIdx) = obj.getScaledPilot();

            residualSignal = rxSignalDaft;
            estimatedPaths = zeros(obj.Config.NumPaths, 3);

            delaySearchRange = 0:obj.Config.PrefixLength;
            dopplerSearchRange = -obj.Config.MaxNormDoppler:obj.Config.MaxNormDoppler;

            % --- SIC 迭代多径估计 ---
            for p = 1:obj.Config.NumPaths
                % 1. 粗搜索
                maxCorr = -inf;
                coarseDelay = 0;
                coarseDoppler = 0;

                for delayTap = delaySearchRange

                    for dopplerShift = dopplerSearchRange
                        [~, cleanRx, energy] = obj.buildAndMaskRx(delayTap, dopplerShift, mockTxFrame);

                        if energy > 1e-10
                            corr = abs(cleanRx' * residualSignal) / sqrt(energy);

                            if corr > maxCorr
                                maxCorr = corr;
                                coarseDelay = delayTap;
                                coarseDoppler = dopplerShift;
                            end

                        end

                    end

                end

                % 精搜索
                costFunc = @(fracDop) -abs(obj.projectGain(coarseDelay, coarseDoppler + fracDop, mockTxFrame, residualSignal));
                [bestFracDoppler, ~] = fminbnd(costFunc, -0.5, 0.5, optimset('TolX', 1e-4));
                finalDoppler = coarseDoppler + bestFracDoppler;

                % 提取增益并执行 SIC
                [~, cleanRx, energy] = obj.buildAndMaskRx(coarseDelay, finalDoppler, mockTxFrame);

                if energy > 1e-10
                    estimatedGain = (cleanRx' * residualSignal) / energy;
                else
                    estimatedGain = 0;
                end

                residualSignal = residualSignal - (estimatedGain * cleanRx);
                estimatedPaths(p, :) = [coarseDelay, finalDoppler, estimatedGain];
            end

            % 全局 LS 优化
            finalBasisMatrix = zeros(obj.Config.NumDataSubcarriers, obj.Config.NumPaths);

            for p = 1:obj.Config.NumPaths
                [~, cleanRx, ~] = obj.buildAndMaskRx(estimatedPaths(p, 1), estimatedPaths(p, 2), mockTxFrame);
                finalBasisMatrix(:, p) = cleanRx;
            end

            refinedGains = pinv(finalBasisMatrix) * rxSignalDaft;
            estimatedPaths(:, 3) = refinedGains;

            % 重建等效信道矩阵
            finalPaths = estimatedPaths;
            physicalChannelMatrix = obj.rebuildLtvChannel(finalPaths);
            effectiveChannelMatrix = AfdmTransforms.generateEffectiveChannelMatrix(physicalChannelMatrix, obj.Config);
        end

        function [expectedRx, maskedRx, energy] = buildAndMaskRx(obj, delayTap, dopplerShift, mockTxFrame)
            tempParams = [delayTap, dopplerShift, 1];
            tempPhysicalChannelMatrix = obj.rebuildLtvChannel(tempParams);
            tempEffectiveChannelMatrix = AfdmTransforms.generateEffectiveChannelMatrix(tempPhysicalChannelMatrix, obj.Config);
            expectedRx = tempEffectiveChannelMatrix * mockTxFrame;

            % 掩膜处理
            numSubcarriers = obj.Config.NumDataSubcarriers;

            if upper(obj.Config.PilotType) == "SINGLE"
                % 单插入式导频以峰值为中心，矩形窗截断
                winRadius = 3;
                [~, peakIndex] = max(abs(expectedRx));
                windowIndicesRaw = (peakIndex - winRadius:peakIndex + winRadius).';
                windowIndices = mod(windowIndicesRaw - 1, numSubcarriers) + 1;
                mask = zeros(numSubcarriers, 1);
                mask(windowIndices) = 1;
                maskedRx = expectedRx .* mask;
            elseif upper(obj.Config.PilotType) == "CAZAC"

                maxAmp = max(abs(expectedRx));
                mask = abs(expectedRx) > (0.1 * maxAmp);
                maskedRx = expectedRx .* mask;
            else
                maskedRx = expectedRx;
            end

            energy = maskedRx' * maskedRx;
        end

        function gain = projectGain(obj, delayTap, dopplerShift, mockTxFrame, rxSignal)
            [~, maskedRx, energy] = obj.buildAndMaskRx(delayTap, dopplerShift, mockTxFrame);

            if energy > 1e-10
                gain = (maskedRx' * rxSignal) / energy;
            else
                gain = 0;
            end

        end

        function physicalChannelMatrix = rebuildLtvChannel(obj, pathParams)
            totalSubcarriers = obj.Config.TotalSubcarriers;
            numDataSubcarriers = obj.Config.NumDataSubcarriers;
            physicalChannelMatrix = sparse(totalSubcarriers, totalSubcarriers);
            numPaths = size(pathParams, 1);

            dopplerCorrectionFactor = totalSubcarriers / numDataSubcarriers;

            for i = 1:numPaths
                pathDelay = round(pathParams(i, 1));
                pathDoppler = pathParams(i, 2) * dopplerCorrectionFactor;
                pathGain = pathParams(i, 3);

                rowIndices = (1 + pathDelay):totalSubcarriers;
                colIndices = 1:(totalSubcarriers - pathDelay);
                delayPermutationMatrix = sparse(rowIndices, colIndices, ones(1, length(rowIndices)), totalSubcarriers, totalSubcarriers);

                timeIndices = (0:totalSubcarriers - 1).';
                dopplerDiagonal = exp(-1j * 2 * pi * pathDoppler * timeIndices / totalSubcarriers);
                dopplerMatrix = spdiags(dopplerDiagonal, 0, totalSubcarriers, totalSubcarriers);

                physicalChannelMatrix = physicalChannelMatrix + pathGain * dopplerMatrix * delayPermutationMatrix;
            end

        end

        function scaledSequence = getScaledPilot(obj)
            baseSeq = obj.Config.PilotSequence;

            if upper(obj.Config.PilotType) == "SINGLE"
                scaledSequence = sqrt(obj.currentPilotPower) * (baseSeq / abs(baseSeq));
            elseif upper(obj.Config.PilotType) == "CAZAC"
                scaledSequence = sqrt(obj.currentPilotPower) * baseSeq;
            else
                scaledSequence = baseSeq;
            end

        end

    end

end
