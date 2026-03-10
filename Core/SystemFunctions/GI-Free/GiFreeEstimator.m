classdef GiFreeEstimator < handle
    % GiFreeEstimator  GI-Free AFDM 信道估计器 (handle 类)
    %
    %   持有 GiFreeConfig 与 FractionalChannelBuilder 的 handle 引用,
    %   所有方法通过 obj 访问共享配置和信道构造能力, 消除参数传递冗余.
    %
    %   整合三大信道估计能力于一体:
    %     (1) 正则化 OMP: 从导频观测中检测路径位置并粗估增益
    %     (2) DD 多普勒精炼: 利用全帧信号精搜分数多普勒偏移
    %     (3) DD 增益重估: 利用全帧信号重估路径复增益
    %   并提供截断均值 (Trimmed Mean) 噪声功率估计工具。
    %
    %   用法:
    %     cfg       = GiFreeConfig();
    %     chBuilder = FractionalChannelBuilder(cfg);
    %     estimator = GiFreeEstimator(cfg, chBuilder);
    %     [estPaths, hEff] = estimator.estimateByOmp(rxSignal);

    properties (SetAccess = private)
        Config % GiFreeConfig handle 引用
        ChannelBuilder % FractionalChannelBuilder handle 引用
    end

    methods

        %% ---------- 构造函数 ----------
        function obj = GiFreeEstimator(cfg, chBuilder)

            arguments
                cfg (1, 1) GiFreeConfig
                chBuilder (1, 1) FractionalChannelBuilder
            end

            obj.Config = cfg;
            obj.ChannelBuilder = chBuilder;
        end

        %% ======================================================================
        %%  正则化 OMP 信道估计
        %% ======================================================================
        function [estPaths, hEffective] = estimateByOmp(obj, rxSignal, ompRegParam)
            % estimateByOmp  干扰感知正则化 OMP 信道估计
            %   rxSignal    : 接收信号 (或经 SIC 清洁后的导频观测)
            %   ompRegParam : 增益估计正则化参数 (可选, 默认 = 0)

            if nargin < 3, ompRegParam = 0; end

            numSc          = obj.Config.NumSubcarriers;
            numPathsUpper  = obj.Config.NumPathsUpper;
            pilotAmplitude = obj.Config.PilotAmplitude;

            estPaths        = zeros(numPathsUpper, 3);
            pilotDictionary = zeros(numSc, numPathsUpper);
            delayRecords    = zeros(numPathsUpper, 1);
            dopplerRecords  = zeros(numPathsUpper, 1);
            residual        = rxSignal;

            for pathIdx = 1:numPathsUpper

                % 步骤 1: 粗搜 — 整数格点峰值检测
                peakHit    = obj.findStrongestPeak(residual);
                delayVal   = peakHit(1);
                intDoppler = peakHit(2);

                % 步骤 2: 精搜 — fminbnd 连续优化分数多普勒
                if obj.Config.SpreadWidth == 0
                    fracDoppler = intDoppler;
                else
                    fracDoppler = obj.refineFractionalDoppler( ...
                        residual, delayVal, intDoppler);
                end

                delayRecords(pathIdx) = delayVal;
                dopplerRecords(pathIdx) = fracDoppler;

                % 步骤 3: 扩展导频字典
                pilotResponse = obj.ChannelBuilder.buildPilotResponseVector( ...
                    delayVal, fracDoppler);
                pilotDictionary(:, pathIdx) = pilotResponse * pilotAmplitude;

                % 步骤 4: 正则化联合 LS 增益估计
                activeDict = pilotDictionary(:, 1:pathIdx);
                gramMatrix = activeDict' * activeDict;

                if ompRegParam > 0
                    jointGains = (gramMatrix + ompRegParam * eye(pathIdx)) \ ...
                        (activeDict' * rxSignal);
                else
                    jointGains = activeDict \ rxSignal;
                end

                % 步骤 5: 正交化残差
                residual = rxSignal - activeDict * jointGains;
            end

            % 组装输出
            for pathIdx = 1:numPathsUpper
                estPaths(pathIdx, :) = [delayRecords(pathIdx), ...
                                            dopplerRecords(pathIdx), ...
                                            jointGains(pathIdx)];
            end

            hEffective = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
        end

        %% ======================================================================
        %%  DD 增益重估 (固定路径位置, 仅更新复增益)
        %% ======================================================================
        function [refinedPaths, hEffective] = refineGainByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)
            % refineGainByDd  利用全帧信号重估路径增益

            if nargin < 5, ddRegParam = 0; end

            numSc    = obj.Config.NumSubcarriers;
            numPaths = size(estPaths, 1);

            if numPaths == 0
                refinedPaths = estPaths;
                hEffective = sparse(numSc, numSc);
                return;
            end

            % 构造全帧响应字典: fullDict(:,i) = H_i * x̂
            fullDict = zeros(numSc, numPaths);

            for i = 1:numPaths
                hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                    estPaths(i, 1), estPaths(i, 2));
                fullDict(:, i) = hSingle * txFrameEst;
            end

            % 正则化 LS 增益重估
            gramFull = fullDict' * fullDict;

            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ (fullDict' * rxSignal);
            else
                newGains = fullDict \ rxSignal;
            end

            refinedPaths = estPaths;
            refinedPaths(:, 3) = newGains;
            hEffective = obj.ChannelBuilder.buildEffectiveChannel(refinedPaths);
        end

        %% ======================================================================
        %%  DD 多普勒精炼 + 增益重估
        %% ======================================================================
        function [refinedPaths, hEffective] = refineDopplerByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)
            % refineDopplerByDd  坐标下降精搜分数多普勒, 然后联合重估增益

            if nargin < 5, ddRegParam = 0; end

            numSc    = obj.Config.NumSubcarriers;
            spreadKv = obj.Config.SpreadWidth;
            numPaths = size(estPaths, 1);

            if numPaths == 0 || spreadKv == 0
                refinedPaths = estPaths;
                hEffective   = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
                return;
            end

            refinedPaths = estPaths;

            % 逐路径多普勒精炼 (坐标下降)
            for i = 1:numPaths
                % 计算排除第 i 条路径后的残差
                residualExcl = rxSignal;

                for j = 1:numPaths
                    if j == i, continue; end
                    hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                        refinedPaths(j, 1), refinedPaths(j, 2));
                    residualExcl = residualExcl - refinedPaths(j, 3) * hSingle * txFrameEst;
                end

                % 全帧匹配滤波搜索最优多普勒
                intDoppler = round(refinedPaths(i, 2));
                delayVal   = refinedPaths(i, 1);

                % 闭包捕获 obj, 使 fminbnd 可调用实例方法
                objFcn = @(nu) -obj.computeFullFrameMetric(delayVal, nu, residualExcl, txFrameEst);

                optOpts            = optimset('TolX', 1e-6, 'Display', 'off');
                bestNu             = fminbnd(objFcn, intDoppler - 0.499, intDoppler + 0.499, optOpts);
                refinedPaths(i, 2) = bestNu;
            end

            % 联合重估所有增益
            fullDict = zeros(numSc, numPaths);

            for i = 1:numPaths
                hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                    refinedPaths(i, 1), refinedPaths(i, 2));
                fullDict(:, i) = hSingle * txFrameEst;
            end

            gramFull = fullDict' * fullDict;

            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ (fullDict' * rxSignal);
            else
                newGains = fullDict \ rxSignal;
            end

            refinedPaths(:, 3) = newGains;

            hEffective = obj.ChannelBuilder.buildEffectiveChannel(refinedPaths);
        end

        %% ======================================================================
        %%  截断均值噪声功率估计
        %% ======================================================================
        function noisePowerEst = estimateNoiseTrimmed(~, rxSignal, hEff, txFrameEst, trimAlpha)
            % estimateNoiseTrimmed  使用截断均值估计等效噪声功率
            %   不依赖 Config, 但保留为实例方法以统一调用风格.
            %   去掉最大的 trimAlpha 比例后取均值, 免疫错误符号的离群效应.

            if nargin < 5, trimAlpha = 0.1; end

            residualPowers = abs(rxSignal - hEff * txFrameEst).^2;
            sortedPowers   = sort(residualPowers, 'ascend');
            numTotal       = length(sortedPowers);

            numKeep       = max(floor(numTotal * (1 - trimAlpha)), 1);
            noisePowerEst = mean(sortedPowers(1:numKeep));
        end

        %% ======================================================================
        %%  内部辅助: 整数格点峰值检测
        %% ======================================================================
        function peakHit = findStrongestPeak(obj, residualSignal)
            numSc          = obj.Config.NumSubcarriers;
            maxDelaySpread = obj.Config.MaxDelaySpread;
            maxDopplerIdx  = obj.Config.MaxDopplerIndex;
            locStep        = obj.Config.LocStep;

            maxMagnitude = -1;
            bestDelay    = 0;
            bestDoppler  = 0;

            for delayIdx = 0:maxDelaySpread

                for dopplerIdx = -maxDopplerIdx:maxDopplerIdx
                    locIndex  = dopplerIdx + locStep * delayIdx;
                    centerBin = mod(-locIndex, numSc);
                    gridVal   = abs(residualSignal(centerBin + 1));

                    if gridVal > maxMagnitude
                        maxMagnitude = gridVal;
                        bestDelay    = delayIdx;
                        bestDoppler  = dopplerIdx;
                    end

                end

            end

            peakHit = [bestDelay, bestDoppler];
        end

        %% ======================================================================
        %%  内部辅助: 分数多普勒连续精搜 (基于导频)
        %% ======================================================================
        function estDoppler = refineFractionalDoppler(obj, residualSignal, delayVal, intDoppler)
            % 闭包捕获 obj, 使 fminbnd 可调用实例方法
            objFcn = @(fracShift) -obj.computePilotMetric( ...
                fracShift, residualSignal, delayVal, intDoppler);

            optOpts      = optimset('TolX', 1e-4, 'Display', 'off');
            fracShiftEst = fminbnd(objFcn, -0.499, 0.499, optOpts);
            estDoppler   = intDoppler + fracShiftEst;
        end

        %% ======================================================================
        %%  内部辅助: 导频匹配滤波度量
        %% ======================================================================
        function metricVal = computePilotMetric(obj, fracShift, signalVec, delayVal, intDoppler)
            pilotResp = obj.ChannelBuilder.buildPilotResponseVector( ...
                delayVal, intDoppler + fracShift);
            normSq = real(pilotResp' * pilotResp);

            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(pilotResp' * signalVec)^2 / normSq;
            end

        end

        %% ======================================================================
        %%  内部辅助: 全帧匹配滤波度量 (DD 多普勒搜索用)
        %% ======================================================================
        function metricVal = computeFullFrameMetric(obj, delayVal, doppler, signalVec, txFrameEst)
            hSingle  = obj.ChannelBuilder.buildPathMatrix(delayVal, doppler);
            response = hSingle * txFrameEst;
            normSq   = real(response' * response);

            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(response' * signalVec)^2 / normSq;
            end

        end

    end

end
