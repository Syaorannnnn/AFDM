classdef FractionalChannelBuilder < handle
    % FractionalChannelBuilder  分数多普勒信道矩阵构造工具 (handle 类)
    %
    %   持有 GiFreeConfig 的引用, 无需在每次调用时传入配置参数.
    %   提供 Dirichlet 系数计算、单径矩阵构造、有效信道合成、
    %   导频响应向量生成，以及随机信道实现生成等方法。
    %
    %   用法:
    %     cfg = GiFreeConfig();
    %     chBuilder = FractionalChannelBuilder(cfg);
    %     hEff = chBuilder.buildEffectiveChannel(pathParams);

    properties (SetAccess = private)
        Config % GiFreeConfig handle 引用
    end

    methods

        %% ---------- 构造函数 ----------
        function obj = FractionalChannelBuilder(cfg)
            % FractionalChannelBuilder  构造函数
            %   cfg : GiFreeConfig handle 对象
            arguments
                cfg (1, 1) GiFreeConfig
            end

            obj.Config = cfg;
        end

        %% ---------- 单径信道矩阵 ----------
        function hPath = buildPathMatrix(obj, delay, doppler)
            % buildPathMatrix  为单条路径 (delay, doppler) 构造 N x N 稀疏信道矩阵
            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            spreadKv = obj.Config.SpreadWidth;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            kvRange  = -spreadKv:spreadKv;
            numKvPts = length(kvRange);

            % 预计算所有 Dirichlet 系数
            dCoeffs = zeros(numKvPts, 1);

            for idx = 1:numKvPts
                dCoeffs(idx) = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvRange(idx), numSc);
            end

            % 稀疏矩阵三元组预分配
            rowIdx = zeros(numSc * numKvPts, 1);
            colIdx = zeros(numSc * numKvPts, 1);
            vals   = complex(zeros(numSc * numKvPts, 1));
            ptr    = 0;

            for p = 0:(numSc - 1)

                for idx = 1:numKvPts
                    if abs(dCoeffs(idx)) < 1e-15, continue; end
                    kvShift     = kvRange(idx);
                    qIdx        = mod(p + locIndex + kvShift, numSc);
                    phaseVal    = exp(1j * 2 * pi / numSc * (numSc * chirpC1 * delay^2 - qIdx * delay + numSc * chirpC2 * (qIdx^2 - p^2)));
                    ptr         = ptr + 1;
                    rowIdx(ptr) = p + 1;
                    colIdx(ptr) = qIdx + 1;
                    vals(ptr)   = phaseVal * dCoeffs(idx);
                end

            end

            hPath = sparse(rowIdx(1:ptr), colIdx(1:ptr), vals(1:ptr), numSc, numSc);
        end

        %% ---------- 多径有效信道 ----------
        function hEff = buildEffectiveChannel(obj, pathParams)
            % buildEffectiveChannel  叠加所有路径，构造总有效信道矩阵
            %   pathParams : P x 3 矩阵, 每行 [delay, doppler, complexGain]
            numSc = obj.Config.NumSubcarriers;
            hEff  = sparse(numSc, numSc);
            if isempty(pathParams), return; end

            for i = 1:size(pathParams, 1)
                hSingle = obj.buildPathMatrix(pathParams(i, 1), pathParams(i, 2));
                hEff    = hEff + pathParams(i, 3) * hSingle;
            end

        end

        %% ---------- 导频响应向量 ----------
        function hCol = buildPilotResponseVector(obj, delay, doppler)
            % buildPilotResponseVector  构造导频子载波 (q=0) 对应的信道列向量
            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            spreadKv = obj.Config.SpreadWidth;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            hCol = zeros(numSc, 1);

            for kvShift = -spreadKv:spreadKv
                pIdx     = mod(-locIndex - kvShift, numSc);
                phaseVal = exp(1j * 2 * pi * (chirpC1 * delay^2 - chirpC2 * pIdx^2));
                dCoeff   = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvShift, numSc);
                hCol(pIdx + 1) = phaseVal * dCoeff;
            end

        end

        %% ---------- 随机信道实现生成 ----------
        function [pathParams, hEff] = generateChannel(obj)
            % generateChannel  按配置随机生成一次信道实现
            cfg           = obj.Config;
            numPaths      = cfg.NumPaths;
            maxDelay      = cfg.MaxDelaySpread;
            maxDopplerIdx = cfg.MaxDopplerIndex;

            % 从 0:maxDelay 中不重复选取时延
            allDelays = 0:maxDelay;
            selectedDelays = allDelays(randperm(length(allDelays), numPaths))';

            % 随机多普勒: 余弦模型 nu_i = k_max * cos(theta_i)
            randomPhases = -pi + 2 * pi * rand(numPaths, 1);
            dopplerShifts = maxDopplerIdx * cos(randomPhases);

            if ~cfg.UseFractionalDoppler
                dopplerShifts = round(dopplerShifts);
            end

            % 确保同一 (delay, round(doppler)) 格点不出现碰撞
            for pathIdx = 2:numPaths
                attempts = 0;

                while any(selectedDelays(1:pathIdx - 1) == selectedDelays(pathIdx) & ...
                        round(dopplerShifts(1:pathIdx - 1)) == round(dopplerShifts(pathIdx)))
                    randomPhases(pathIdx)  = -pi + 2 * pi * rand();
                    dopplerShifts(pathIdx) = maxDopplerIdx * cos(randomPhases(pathIdx));

                    if ~cfg.UseFractionalDoppler
                        dopplerShifts(pathIdx) = round(dopplerShifts(pathIdx));
                    end

                    attempts = attempts + 1;
                    if attempts > 100, break; end
                end

            end

            % 等功率 Rayleigh 增益
            complexGains = sqrt(1 / (2 * numPaths)) * ...
                (randn(numPaths, 1) + 1j * randn(numPaths, 1));
            pathParams = [selectedDelays, dopplerShifts, complexGains];
            hEff       = obj.buildEffectiveChannel(pathParams);
        end

    end

    %% ========== 纯数学工具 — 保留为 Static ==========
    methods (Static)

        function dVal = computeDirichletCoeff(fracPart, shiftIndex, numSc)
            % computeDirichletCoeff  计算分数多普勒的 Dirichlet 核系数
            %   此函数不依赖任何系统配置, 属于纯数学运算, 故保留 Static.
            xVal = fracPart - shiftIndex;

            if abs(xVal) < 1e-12
                dVal = 1;
            elseif abs(fracPart) < 1e-12
                dVal = 0;
            else
                dVal = (1 / numSc) * (1 - exp(-1j * 2 * pi * fracPart)) / ...
                    (1 - exp(-1j * 2 * pi * xVal / numSc));
            end

        end

    end

end
