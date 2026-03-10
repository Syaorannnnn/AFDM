classdef EpReceiver < handle
    % EpReceiver: EP-AFDM 接收机 (单导频版, 修复版)
    %
    % 修复记录:
    %   Bug#1 — rebuildChannel 现在构造循环矩阵, 与 LtvChannel 结构匹配.
    %           原版只构建了三角矩阵, 缺失环绕项, 导致等效信道完全错误.
    %   Bug#2 — IP2D 消除现在始终执行 (不再限于 CAZAC 模式).
    %           导频功率远大于数据, 不消除则数据被导频干扰淹没.
    %
    % 接收流程 (单导频 + ZP 保护带):
    %   1. 去 CPP → DAFT 变换
    %   2. 信道估计 (时域相关搜 + 分数多普勒精搜 + 全局 LS)
    %   3. IP2D 消除: 从接收信号中减去导频的信道响应
    %   4. MMSE 均衡数据子载波
    %   5. 硬判决
    %
    % 由于 ZP 保护带保证了导频观测不受数据干扰,
    % 信道估计本身就是干净的, 无需迭代 SIC.

    properties (Access = public)
        CsiMode      (1, 1) string = "Estimated"  % "Perfect" 或 "Estimated"
        EqualizerType (1, 1) string = "MMSE"
    end

    properties (Access = private)
        Config
        currentPilotPower (1, 1) double = 1
    end

    methods (Access = public)

        function obj = EpReceiver(configObj)
            obj.Config = configObj;
        end

        function rxData = receive(obj, rxSignal, noisePower, physicalChannelMatrix, pilotPower)
            obj.currentPilotPower = pilotPower;

            % 去 CPP 前缀
            rxNoPrefix = rxSignal(obj.Config.PrefixLength + 1:end);

            % DAFT 变换
            if upper(obj.Config.WaveformType) == "AFDM"
                daftRx = AfdmTransforms.daft(rxNoPrefix, ...
                    obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                daftRx = AfdmTransforms.dft(rxNoPrefix);
            end

            % 获取等效信道矩阵
            if obj.CsiMode == "Perfect"
                hEff = AfdmTransforms.generateEffectiveChannelMatrix( ...
                    physicalChannelMatrix, obj.Config);
            else
                [hEff, ~] = obj.runChannelEstimator(daftRx);
            end

            % ======================================================
            %  IP2D 消除 (Bug#2 修复: 始终执行)
            %
            %  导频功率 >> 数据功率, 若不消除, 导频通过信道扩散到
            %  所有子载波上的能量会完全淹没数据信号.
            %  做法: 构造仅含导频的 DAFT 域帧, 通过等效信道矩阵
            %        计算其对所有子载波的贡献, 从接收信号中减去.
            % ======================================================
            pilotFrame = zeros(obj.Config.NumDataSubcarriers, 1);
            pilotFrame(obj.Config.PilotIndex) = sqrt(obj.currentPilotPower);
            pilotContribution = hEff * pilotFrame;
            cleanDataSignal = daftRx - pilotContribution;

            % MMSE 均衡 (仅作用于数据子载波列)
            hData = hEff(:, obj.Config.ActiveIndices);
            estData = (hData' * hData + noisePower * eye(obj.Config.NumActiveCarriers)) ...
                \ (hData' * cleanDataSignal);

            % 硬判决
            rxData = qamdemod(estData, obj.Config.ModulationOrder, 'UnitAveragePower', true);
            rxData = rxData(:);
        end

    end

    methods (Access = private)

        % ================================================================
        %  信道估计器: 时域相关搜索 + 分数多普勒精搜 + 全局 LS
        % ================================================================
        function [hEff, finalPaths] = runChannelEstimator(obj, rxSignalDaft)
            N = obj.Config.NumDataSubcarriers;
            c1 = obj.Config.ChirpParam1;
            c2 = obj.Config.ChirpParam2;

            % 构造仅含导频的参考帧
            mockTxFrame = zeros(N, 1);
            mockTxFrame(obj.Config.PilotIndex) = sqrt(obj.currentPilotPower);

            % 预计算时域导频
            if upper(obj.Config.WaveformType) == "AFDM"
                timePilot = AfdmTransforms.idaft(mockTxFrame, c1, c2);
            else
                timePilot = AfdmTransforms.idft(mockTxFrame);
            end

            residual = rxSignalDaft;
            estPaths = zeros(obj.Config.NumPaths, 3);

            delayRange = 0:obj.Config.MaxPathDelays;
            dopplerRange = -obj.Config.MaxNormDoppler:obj.Config.MaxNormDoppler;

            % 预计算多普勒相位矩阵 (向量化粗搜)
            timeVec = (0:N - 1).';
            normDopplerVec = dopplerRange(:).' / N;
            dopplerPhaseMat = exp(-1j * 2 * pi * timeVec .* normDopplerVec);

            for pathIdx = 1:obj.Config.NumPaths
                bestCorr = -inf;
                bestDelay = 0;
                bestDopIdx = 1;

                % --- 向量化粗搜: 遍历时延, 批量处理多普勒 ---
                for dIdx = 1:length(delayRange)
                    curDelay = delayRange(dIdx);
                    shifted = circshift(timePilot, curDelay);
                    candidates = shifted .* dopplerPhaseMat;  % N × numDoppler

                    if upper(obj.Config.WaveformType) == "AFDM"
                        daftCand = AfdmTransforms.daft(candidates, c1, c2);
                    else
                        daftCand = AfdmTransforms.dft(candidates);
                    end

                    % 窗口掩膜: 仅保留峰值附近 ±3 bin (单导频特征)
                    for col = 1:size(daftCand, 2)
                        daftCand(:, col) = obj.applyPeakWindow(daftCand(:, col), 3);
                    end

                    energies = sum(abs(daftCand) .^ 2, 1);
                    corrs = abs(residual' * daftCand) ./ sqrt(energies + 1e-12);
                    [maxVal, maxIdx] = max(corrs);

                    if maxVal > bestCorr
                        bestCorr = maxVal;
                        bestDelay = curDelay;
                        bestDopIdx = maxIdx;
                    end
                end

                coarseDoppler = dopplerRange(bestDopIdx);

                % --- 分数多普勒精搜 ---
                costFn = @(frac) -obj.calcCorrelation( ...
                    bestDelay, coarseDoppler + frac, timePilot, residual, [c1, c2]);
                [fracOpt, ~] = fminbnd(costFn, -0.5, 0.5, optimset('TolX', 1e-3));
                finalDoppler = coarseDoppler + fracOpt;

                % --- 提取当前径, 更新残差 (SIC) ---
                [cleanSig, pathEnergy] = obj.simulatePathEffect( ...
                    bestDelay, finalDoppler, timePilot, [c1, c2]);
                estGain = (cleanSig' * residual) / (pathEnergy + 1e-12);
                residual = residual - estGain * cleanSig;
                estPaths(pathIdx, :) = [bestDelay, finalDoppler, estGain];
            end

            % --- 全局 LS 联合优化增益 ---
            basisMat = zeros(N, obj.Config.NumPaths);

            for pathIdx = 1:obj.Config.NumPaths
                basisMat(:, pathIdx) = obj.simulatePathEffect( ...
                    estPaths(pathIdx, 1), estPaths(pathIdx, 2), timePilot, [c1, c2]);
            end

            refinedGains = lsqminnorm(basisMat, rxSignalDaft);
            estPaths(:, 3) = refinedGains;
            finalPaths = estPaths;

            % --- 重建等效信道 ---
            physChannel = obj.rebuildChannel(finalPaths);
            hEff = AfdmTransforms.generateEffectiveChannelMatrix(physChannel, obj.Config);
        end

        % ================================================================
        %  模拟单径对 DAFT 域导频的效果 (用于相关搜索)
        % ================================================================
        function [daftSig, energy] = simulatePathEffect(obj, delay, doppler, timePilot, chirpParams)
            N = size(timePilot, 1);
            timeVec = (0:N - 1).';
            phaseRot = exp(-1j * 2 * pi * doppler / N .* timeVec);

            pathTimeDomain = circshift(timePilot, delay) .* phaseRot;

            if upper(obj.Config.WaveformType) == "AFDM"
                daftSig = AfdmTransforms.daft(pathTimeDomain, chirpParams(1), chirpParams(2));
            else
                daftSig = AfdmTransforms.dft(pathTimeDomain);
            end

            % 峰值窗口掩膜
            daftSig = obj.applyPeakWindow(daftSig, 3);

            if nargout > 1
                energy = sum(abs(daftSig) .^ 2);
            end
        end

        % ================================================================
        %  峰值窗口: 保留峰值附近 ±radius 的 bin, 其余置零
        % ================================================================
        function masked = applyPeakWindow(~, sig, radius)
            N = length(sig);
            [~, peakIdx] = max(abs(sig));
            winRaw = (peakIdx - radius:peakIdx + radius).';
            winIdx = mod(winRaw - 1, N) + 1;
            mask = zeros(N, 1);
            mask(winIdx) = 1;
            masked = sig .* mask;
        end

        % ================================================================
        %  相关性代价函数 (分数多普勒精搜用)
        % ================================================================
        function val = calcCorrelation(obj, delay, doppler, timePilot, rxSig, chirpParams)
            [daftSig, energy] = obj.simulatePathEffect(delay, doppler, timePilot, chirpParams);

            if energy > 1e-12
                val = abs(daftSig' * rxSig) / sqrt(energy);
            else
                val = 0;
            end
        end

        % ================================================================
        %  从路径参数重建物理信道矩阵 (Bug#1 修复: 循环结构)
        %
        %  LtvChannel 构造的是循环移位 + 多普勒相位矩阵:
        %    H(r, c) = Σ_i  h_i · exp(-j2π·ν_i·r/N_tot)
        %              当 r = mod(c + l_i, N_tot) 时
        %
        %  原版 rebuildChannel 只生成了 r > l_i 的下三角部分,
        %  缺少了 r < l_i 的环绕部分 (当 c + l_i ≥ N_tot 时 r 应回绕到 0).
        %  这导致重建矩阵与真实信道在结构上完全不同.
        %
        %  修复: 用 mod 运算正确生成循环移位的行索引.
        % ================================================================
        function physChannel = rebuildChannel(obj, pathParams)
            totalSc = obj.Config.TotalSubcarriers;
            numDataSc = obj.Config.NumDataSubcarriers;
            numPaths = size(pathParams, 1);

            % 多普勒缩放: 估计器用 N_data 归一化, 物理信道用 N_total 归一化
            dopplerScale = totalSc / numDataSc;

            allRows = [];
            allCols = [];
            allVals = [];

            colVec = (0:totalSc - 1).';  % 0-based 列索引

            for i = 1:numPaths
                delay = round(pathParams(i, 1));
                doppler = pathParams(i, 2) * dopplerScale;
                gain = pathParams(i, 3);

                % 循环移位: 列 c → 行 mod(c + delay, N_total)
                rowVec = mod(colVec + delay, totalSc);

                % 多普勒相位: 作用在行索引上
                phases = exp(-1j * 2 * pi * doppler * rowVec / totalSc);
                vals = gain * phases;

                % 转为 1-based 索引
                allRows = [allRows; rowVec + 1]; %#ok<AGROW>
                allCols = [allCols; colVec + 1]; %#ok<AGROW>
                allVals = [allVals; vals];        %#ok<AGROW>
            end

            physChannel = sparse(allRows, allCols, allVals, totalSc, totalSc);
        end

    end

end