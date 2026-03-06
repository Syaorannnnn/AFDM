classdef AfdmTransmitter < handle
    % AfdmTransmitter: 负责 AFDM 信号的生成、调制、导频插入与组帧

    properties (Access = private)
        Config % 保存对全局配置对象的引用
    end

    methods (Access = public)

        % --- 构造函数 ---
        % 依赖注入 Config 对象
        function obj = AfdmTransmitter(configObj)
            obj.Config = configObj;
        end

        % --- 发送主函数 ---
        function [txSignal, originalData] = transmit(obj, targetPilotPower)
            % 生成随机原始数据并进行 QAM 调制
            originalData = randi([0, obj.Config.ModulationOrder - 1], obj.Config.NumActiveCarriers, 1);
            qamSymbols = qammod(originalData, obj.Config.ModulationOrder, 'UnitAveragePower', true);

            % 构建 DAFT/DFT 域的帧
            freqFrame = zeros(obj.Config.NumDataSubcarriers, 1);

            % 插入并缩放导频 (根据不同导频策略处理功率)
            scaledPilot = obj.scalePilot(targetPilotPower);
            startIdx = obj.Config.PilotIndex;
            endIdx = startIdx + obj.Config.PilotSequenceLength - 1;
            freqFrame(startIdx:endIdx) = scaledPilot;

            % 映射有效数据
            freqFrame(obj.Config.ActiveIndices) = qamSymbols;

            % 脉冲成形
            shapedFreqFrame = freqFrame .* obj.Config.PulseShapingWindow;

            % 时域变换
            if upper(obj.Config.WaveformType) == "AFDM"
                timeFrame = AfdmTransforms.idaft(shapedFreqFrame, obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                timeFrame = AfdmTransforms.idft(freqFrame);
            end

            % 添加前缀 (CP 或 CPP)
            txSignal = obj.addPrefix(timeFrame);
        end

    end

    methods (Access = private)

        % --- 导频功率缩放 ---
        function scaledSequence = scalePilot(obj, targetPower)
            baseSeq = obj.Config.PilotSequence;

            if upper(obj.Config.PilotType) == "SINGLE"
                % 单符号：保证其模长等于 sqrt(targetPower)
                scaledSequence = sqrt(targetPower) * (baseSeq / abs(baseSeq));
            elseif upper(obj.Config.PilotType) == "CAZAC"
                % CAZAC序列：直接乘上幅度因子
                scaledSequence = sqrt(targetPower) * baseSeq;
            else
                scaledSequence = baseSeq;
            end

        end

        % --- 添加前缀 ---
        function txSignal = addPrefix(obj, timeFrame)
            prefixLength = obj.Config.PrefixLength;
            numDataSubcarriers = obj.Config.NumDataSubcarriers;

            if upper(obj.Config.WaveformType) == "AFDM"
                % AFDM 需要添加 CPP (Chirp Periodic Prefix)
                phaseVector = exp(-1j * 2 * pi * obj.Config.ChirpParam1 * (numDataSubcarriers ^ 2 + 2 * numDataSubcarriers * (-prefixLength:-1).'));
                preFix = timeFrame(end - prefixLength + 1:end) .* phaseVector;
            else
                % 传统 OFDM 的 CP
                preFix = timeFrame(end - prefixLength + 1:end);
            end

            txSignal = [preFix; timeFrame];
        end

    end

end
