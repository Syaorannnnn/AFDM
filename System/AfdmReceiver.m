classdef AfdmReceiver < handle
    % AfdmReceiver: AFDM 接收机 — 整理版
    %
    % 核心功能:
    %   1. 信道估计: 时域相关搜索 + 分数多普勒精搜 + 全局 LS 优化
    %   2. 均衡器: MMSE (线性) / MRC-DFE (迭代非线性)
    %   3. 迭代 SIC: CAZAC 导频下的导频-数据互干扰消除
    %
    % 迭代策略 (自适应收敛):
    %   外层 SIC 和内层 MRC-DFE 均采用 "最大迭代上限 + 收敛阈值":
    %     - 每轮检查残差/符号变化量, 低于阈值则提前退出
    %     - 高 SNR: 2~3 轮即收敛, 节省计算
    %     - 低 SNR: 可能跑满上限, 保证充分收敛
    %
    % 已移除: 脉冲成形窗 / 匹配滤波 / WMRC 均衡器 (后续可扩展)

    properties (Access = public)
        CsiMode (1, 1) string = "Estimated" % "Perfect" 或 "Estimated"
        EqualizerType (1, 1) string = "MMSE" % "MMSE" 或 "MRC"

        % --- 外层 SIC 配置 ---
        EnableIterativeIc (1, 1) logical = true
        MaxSicIterations (1, 1) double = 8 % 最大轮次上限
        SicConvergenceThreshold (1, 1) double = 1e-4 % 归一化符号变化量阈值

        % --- 内层 MRC-DFE 配置 ---
        MaxMrcIterations (1, 1) double = 20 % 最大轮次上限
        MrcConvergenceThreshold (1, 1) double = 1e-6 % 估计向量范数变化阈值
    end

    properties (Access = private)
        Config
        currentPilotPower (1, 1) double = 1;
    end

    methods (Access = public)

        function obj = AfdmReceiver(configObj)
            obj.Config = configObj;
        end

        function rxData = receive(obj, rxSignal, noisePower, physicalChannelMatrix, pilotPower)
            obj.currentPilotPower = pilotPower;

            % 去前缀
            rxNoPrefix = rxSignal(obj.Config.PrefixLength + 1:end);

            % DAFT / DFT 变换
            if upper(obj.Config.WaveformType) == "AFDM"
                daftRx = AfdmTransforms.daft(rxNoPrefix, obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                daftRx = AfdmTransforms.dft(rxNoPrefix);
            end

            rxData = obj.runDetection(daftRx, noisePower, physicalChannelMatrix);
        end

    end

    methods (Access = private)

        % ================================================================
        % 检测主循环: 自适应 SIC 迭代
        % ================================================================
        function rxData = runDetection(obj, daftRx, noisePower, physicalChannelMatrix)
            residualSignalForPilot = daftRx;
            qamConstellation = qammod(0:obj.Config.ModulationOrder - 1, obj.Config.ModulationOrder, 'UnitAveragePower', true);
            maxAmp = max(abs(qamConstellation));

            isSicEnabled = (obj.CsiMode == "Estimated") && ...
                (upper(obj.Config.PilotType) == "CAZAC") && obj.EnableIterativeIc;

            if isSicEnabled
                maxIters = obj.MaxSicIterations;
            else
                maxIters = 1; % 非 SIC 模式只跑一次
            end

            prevEqSignal = [];
            eqSignal = [];

            for iter = 1:maxIters

                % --- 信道获取 ---
                if obj.CsiMode == "Perfect"
                    effectiveChannelMatrix = AfdmTransforms.generateEffectiveChannelMatrix(physicalChannelMatrix, obj.Config);
                else
                    [effectiveChannelMatrix, ~] = obj.runChannelEstimator(residualSignalForPilot);
                end

                % --- 消除导频对数据的干扰 (IP2D) ---
                if isSicEnabled
                    mockPilotFrame = zeros(obj.Config.NumDataSubcarriers, 1);
                    startIdx = obj.Config.PilotIndex;
                    endIdx = startIdx + obj.Config.PilotSequenceLength - 1;
                    mockPilotFrame(startIdx:endIdx) = obj.getScaledPilot();
                    ip2dInterference = effectiveChannelMatrix * mockPilotFrame;
                    residualSignalForData = daftRx - ip2dInterference;
                else
                    residualSignalForData = daftRx;
                end

                % --- 数据均衡 ---
                eqSignal = obj.runEqualizer(residualSignalForData, effectiveChannelMatrix, noisePower);

                % --- SIC 收敛检测 ---
                % 归一化变化量 = ||当前输出 - 上一轮输出|| / ||当前输出||
                % 低于阈值说明新一轮几乎没改善, 可以提前退出
                if ~isempty(prevEqSignal)
                    sicDelta = norm(eqSignal - prevEqSignal) / max(norm(eqSignal), 1e-12);

                    if sicDelta < obj.SicConvergenceThreshold
                        break;
                    end

                end

                prevEqSignal = eqSignal;

                % --- 迭代重构: 消除数据对导频的干扰 (ID2P) ---
                if iter < maxIters && isSicEnabled

                    if iter == 1
                        % 第一轮: 软判决 + 限幅 (防止错误传播)
                        alpha = 0.6;
                        reModData = eqSignal;
                        exceedIdx = abs(reModData) > maxAmp;
                        reModData(exceedIdx) = maxAmp * exp(1j * angle(reModData(exceedIdx)));
                    else
                        % 后续轮: 硬判决重构 (估计更可靠)
                        alpha = 1.0;
                        tempDemod = qamdemod(eqSignal, obj.Config.ModulationOrder, 'UnitAveragePower', true);
                        reModData = qammod(tempDemod, obj.Config.ModulationOrder, 'UnitAveragePower', true);
                    end

                    mockDataFrame = zeros(obj.Config.NumDataSubcarriers, 1);
                    mockDataFrame(obj.Config.ActiveIndices) = reModData(:);
                    id2pInterference = effectiveChannelMatrix * mockDataFrame;
                    residualSignalForPilot = daftRx - alpha * id2pInterference;
                end

            end

            % 最终硬判决
            rxData = qamdemod(eqSignal, obj.Config.ModulationOrder, 'UnitAveragePower', true);
            rxData = rxData(:);
        end

        % ================================================================
        % 均衡器分发
        % ================================================================
        function eqSignal = runEqualizer(obj, rxSignal, effectiveChannelMatrix, noisePower)

            switch upper(obj.EqualizerType)
                case "MMSE"
                    eqSignal = obj.equalizeMMSE(rxSignal, effectiveChannelMatrix, noisePower);
                case "MRC"
                    eqSignal = obj.equalizeMRC(rxSignal, effectiveChannelMatrix, noisePower);
                otherwise
                    error('AfdmReceiver:UnknownEqualizer', '未知均衡器: %s', obj.EqualizerType);
            end

        end

        % ================================================================
        % MMSE 线性均衡器
        % ================================================================
        function estData = equalizeMMSE(obj, rxSignal, effectiveChannelMatrix, noisePower)
            activeChannelMatrix = effectiveChannelMatrix(:, obj.Config.ActiveIndices);
            gramMatrix = activeChannelMatrix' * activeChannelMatrix;
            estData = (gramMatrix + noisePower * eye(obj.Config.NumActiveCarriers)) \ (activeChannelMatrix' * rxSignal);
        end

        % ================================================================
        % MRC-DFE 迭代均衡器 (自适应收敛)
        %
        % 工作原理:
        %   每轮按信号能量从强到弱逐符号处理:
        %     1. MRC 合并: 从残差中提取当前符号的最佳估计
        %     2. 干扰消除: 用估计值更新残差, 减少对后续符号的干扰
        %   重复直到收敛或达到最大迭代次数
        %
        % 优化点 (相对原版):
        %   - 收敛检测从 norm(全向量差) 改为在内循环中累加增量平方和,
        %     避免每轮额外遍历全向量
        %   - 稀疏矩阵结构预提取到 cell 数组, 消除循环内 find() 开销
        % ================================================================
        function estData = equalizeMRC(obj, rxSignal, effectiveChannelMatrix, noisePower)
            numSubcarriers = obj.Config.NumDataSubcarriers;
            columnEnergies = full(sum(abs(effectiveChannelMatrix) .^ 2, 1)).';

            % 预提取稀疏结构 (每列的非零行索引和值)
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

            % 按列能量从强到弱排序活跃子载波 (强信号先检测, 提高消除质量)
            activeEnergies = columnEnergies(obj.Config.ActiveIndices);
            [~, sortIdx] = sort(activeEnergies, 'descend');
            orderedIndices = obj.Config.ActiveIndices(sortIdx);

            for n = 1:obj.MaxMrcIterations

                % 累加增量平方和, 用于收敛检测 (避免额外的 norm 计算)
                totalDeltaSq = 0;

                for k = orderedIndices
                    rows = columnRowIndices{k};
                    if isempty(rows), continue; end
                    channelValues = columnValues{k};

                    % MRC 合并: 从残差中提取符号 k 的估计
                    mrcOutput = sum(conj(channelValues) .* residualSignal(rows)) ...
                        + columnEnergies(k) * previousEstimate(k);
                    newEstimate = mrcOutput / (columnEnergies(k) + noisePower);

                    % 更新残差 (消除符号 k 的干扰贡献)
                    estimateChange = newEstimate - previousEstimate(k);
                    currentEstimate(k) = newEstimate;

                    if abs(estimateChange) > 1e-6
                        residualSignal(rows) = residualSignal(rows) - channelValues * estimateChange;
                    end

                    totalDeltaSq = totalDeltaSq + abs(estimateChange) ^ 2;
                end

                % 自适应收敛检测 (基于累加增量, 比 norm 更高效)
                if sqrt(totalDeltaSq) < obj.MrcConvergenceThreshold
                    break;
                end

                previousEstimate = currentEstimate;
            end

            estData = currentEstimate(obj.Config.ActiveIndices);
        end

        % ================================================================
        % 信道估计器: 时域相关搜索 + 分数多普勒精搜 + 全局 LS
        % ================================================================
        function [effectiveChannelMatrix, finalPaths] = runChannelEstimator(obj, rxSignalDaft)
            N = obj.Config.NumDataSubcarriers;
            c1 = obj.Config.ChirpParam1;
            c2 = obj.Config.ChirpParam2;

            % 构造仅含导频的参考帧
            mockTxFrame = zeros(N, 1);
            startIndex = obj.Config.PilotIndex;
            endIndex = startIndex + obj.Config.PilotSequenceLength - 1;
            mockTxFrame(startIndex:endIndex) = obj.getScaledPilot();

            % 预计算时域导频
            if upper(obj.Config.WaveformType) == "AFDM"
                timeDomainPilotBase = AfdmTransforms.idaft(mockTxFrame, c1, c2);
            else
                timeDomainPilotBase = AfdmTransforms.idft(mockTxFrame);
            end

            residualSignal = rxSignalDaft;
            estimatedPaths = zeros(obj.Config.NumPaths, 3);

            delaySearchRange = 0:obj.Config.PrefixLength;
            dopplerSearchRange = -obj.Config.MaxNormDoppler:obj.Config.MaxNormDoppler;

            % 预计算多普勒相位矩阵
            timeVector = (0:N - 1).';
            normalizedDopplerVector = dopplerSearchRange(:).' / N;
            dopplerMatrix = exp(-1j * 2 * pi * timeVector .* normalizedDopplerVector);

            for pathIndex = 1:obj.Config.NumPaths
                maxCorrelation = -inf;
                bestDelay = 0;
                bestDopplerIndex = 1;

                % 向量化粗搜
                for delayIdx = 1:length(delaySearchRange)
                    currentDelay = delaySearchRange(delayIdx);

                    % 时域移位 + 批量施加多普勒
                    timeDomainShifted = circshift(timeDomainPilotBase, currentDelay);
                    timeDomainCandidates = timeDomainShifted .* dopplerMatrix;

                    % 批量 DAFT 变换
                    if upper(obj.Config.WaveformType) == "AFDM"
                        daftCandidates = AfdmTransforms.daft(timeDomainCandidates, c1, c2);
                    else
                        daftCandidates = AfdmTransforms.dft(timeDomainCandidates);
                    end

                    % CAZAC 掩膜: 截取幅度 > 10% 峰值的部分
                    if upper(obj.Config.PilotType) == "CAZAC"
                        maxAmplitudes = max(abs(daftCandidates), [], 1);
                        validMask = abs(daftCandidates) > (0.1 * maxAmplitudes);
                        daftCandidates = daftCandidates .* validMask;
                    end

                    % 批量相关性计算
                    candidateEnergies = sum(abs(daftCandidates) .^ 2, 1);
                    correlations = abs(residualSignal' * daftCandidates) ./ sqrt(candidateEnergies +1e-12);

                    [currentMaxVal, currentMaxIdx] = max(correlations);

                    if currentMaxVal > maxCorrelation
                        maxCorrelation = currentMaxVal;
                        bestDelay = currentDelay;
                        bestDopplerIndex = currentMaxIdx;
                    end

                end

                coarseDoppler = dopplerSearchRange(bestDopplerIndex);
                chirpParams = [c1, c2];

                % 分数多普勒精搜
                costFunction = @(fracDoppler) -obj.calcCorrelation( ...
                    bestDelay, coarseDoppler + fracDoppler, timeDomainPilotBase, residualSignal, chirpParams);
                [bestFractionalDoppler, ~] = fminbnd(costFunction, -0.5, 0.5, optimset('TolX', 1e-3));
                finalDoppler = coarseDoppler + bestFractionalDoppler;

                % 提取当前径并 SIC
                [cleanReceivedSignal, pathEnergy] = obj.applyChannelEffect( ...
                    bestDelay, finalDoppler, timeDomainPilotBase, chirpParams);
                estimatedGain = (cleanReceivedSignal' * residualSignal) / (pathEnergy +1e-12);

                residualSignal = residualSignal - (estimatedGain * cleanReceivedSignal);
                estimatedPaths(pathIndex, :) = [bestDelay, finalDoppler, estimatedGain];
            end

            % 全局 LS 优化增益
            finalBasisMatrix = zeros(N, obj.Config.NumPaths);
            chirpParams = [c1, c2];

            for pathIndex = 1:obj.Config.NumPaths
                finalBasisMatrix(:, pathIndex) = obj.applyChannelEffect( ...
                    estimatedPaths(pathIndex, 1), estimatedPaths(pathIndex, 2), timeDomainPilotBase, chirpParams);
            end

            refinedGains = lsqminnorm(finalBasisMatrix, rxSignalDaft);
            estimatedPaths(:, 3) = refinedGains;
            finalPaths = estimatedPaths;

            % 重建等效信道矩阵
            physicalChannelMatrix = obj.rebuildChannel(finalPaths);
            effectiveChannelMatrix = AfdmTransforms.generateEffectiveChannelMatrix(physicalChannelMatrix, obj.Config);
        end

        % --- 模拟单径信道效果 ---
        function [cleanReceivedDaft, signalEnergy] = applyChannelEffect(obj, delayTap, dopplerShift, timeDomainPilot, chirpParams)
            signalLength = size(timeDomainPilot, 1);
            timeVector = (0:signalLength - 1).';
            dopplerRadian = 2 * pi * dopplerShift / signalLength;
            phaseRotationVector = exp(-1j * dopplerRadian .* timeVector);

            pathSignalTimeDomain = circshift(timeDomainPilot, delayTap) .* phaseRotationVector;

            if upper(obj.Config.WaveformType) == "AFDM"
                cleanReceivedDaft = AfdmTransforms.daft(pathSignalTimeDomain, chirpParams(1), chirpParams(2));
            else
                cleanReceivedDaft = AfdmTransforms.dft(pathSignalTimeDomain);
            end

            % 掩膜处理
            if upper(obj.Config.PilotType) == "SINGLE"
                windowRadius = 3;
                [~, peakIndex] = max(abs(cleanReceivedDaft));
                windowIndicesRaw = (peakIndex - windowRadius:peakIndex + windowRadius).';
                windowIndices = mod(windowIndicesRaw - 1, signalLength) + 1;
                validMask = zeros(signalLength, 1);
                validMask(windowIndices) = 1;
                cleanReceivedDaft = cleanReceivedDaft .* validMask;
            elseif upper(obj.Config.PilotType) == "CAZAC"
                maxAmplitude = max(abs(cleanReceivedDaft), [], 1);
                validMask = abs(cleanReceivedDaft) > (0.1 * maxAmplitude);
                cleanReceivedDaft = cleanReceivedDaft .* validMask;
            end

            if nargout > 1
                signalEnergy = sum(abs(cleanReceivedDaft) .^ 2, 1);
            end

        end

        % --- 相关性计算 (精搜代价函数) ---
        function val = calcCorrelation(obj, delayTap, dopplerShift, timeDomainPilot, receivedSignal, chirpParams)
            [cleanReceivedDaft, signalEnergy] = obj.applyChannelEffect(delayTap, dopplerShift, timeDomainPilot, chirpParams);

            if signalEnergy > 1e-12
                val = abs(cleanReceivedDaft' * receivedSignal) / sqrt(signalEnergy);
            else
                val = 0;
            end

        end

        % --- 从路径参数重建 LTV 物理信道矩阵 ---
        function physicalChannelMatrix = rebuildChannel(obj, pathParams)
            totalSubcarriers = obj.Config.TotalSubcarriers;
            numDataSubcarriers = obj.Config.NumDataSubcarriers;
            numPaths = size(pathParams, 1);
            dopplerCorrectionFactor = totalSubcarriers / numDataSubcarriers;

            % 预分配所有路径的行列索引和值, 最后一次性构造稀疏矩阵
            allRows = [];
            allCols = [];
            allVals = [];

            for i = 1:numPaths
                pathDelay = round(pathParams(i, 1));
                pathDoppler = pathParams(i, 2) * dopplerCorrectionFactor;
                pathGain = pathParams(i, 3);

                rowIndices = ((1 + pathDelay):totalSubcarriers).';
                colIndices = (1:(totalSubcarriers - pathDelay)).';

                % 对应行的多普勒相位
                dopplerPhases = exp(-1j * 2 * pi * pathDoppler * (rowIndices - 1) / totalSubcarriers);
                vals = pathGain * dopplerPhases;

                allRows = [allRows; rowIndices]; %#ok<AGROW>
                allCols = [allCols; colIndices]; %#ok<AGROW>
                allVals = [allVals; vals]; %#ok<AGROW>
            end

            physicalChannelMatrix = sparse(allRows, allCols, allVals, totalSubcarriers, totalSubcarriers);
        end

        % --- 导频功率缩放 ---
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
