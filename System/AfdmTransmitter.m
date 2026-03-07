classdef AfdmTransmitter < handle
    % AfdmTransmitter: AFDM 信号生成、调制、导频插入与组帧
    %
    % 整理版 — 清除脉冲成形相关代码, 保留核心发射链路:
    %   QAM 调制 → 导频插入 → 帧组装 → IDAFT 变换 → 添加前缀 (CPP)

    properties (Access = private)
        Config
    end

    methods (Access = public)

        function obj = AfdmTransmitter(configObj)
            obj.Config = configObj;
        end

        function [txSignal, originalData] = transmit(obj, targetPilotPower)
            % 生成随机数据并 QAM 调制
            originalData = randi([0, obj.Config.ModulationOrder - 1], obj.Config.NumActiveCarriers, 1);
            qamSymbols = qammod(originalData, obj.Config.ModulationOrder, 'UnitAveragePower', true);

            % 组装 DAFT 域帧
            freqFrame = zeros(obj.Config.NumDataSubcarriers, 1);

            % 插入导频
            scaledPilot = obj.scalePilot(targetPilotPower);
            startIdx = obj.Config.PilotIndex;
            endIdx = startIdx + obj.Config.PilotSequenceLength - 1;
            freqFrame(startIdx:endIdx) = scaledPilot;

            % 映射数据
            freqFrame(obj.Config.ActiveIndices) = qamSymbols;

            % IDAFT / IDFT 变换到时域
            if upper(obj.Config.WaveformType) == "AFDM"
                timeFrame = AfdmTransforms.idaft(freqFrame, obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                timeFrame = AfdmTransforms.idft(freqFrame);
            end

            % 添加前缀
            txSignal = obj.addPrefix(timeFrame);
        end

    end

    methods (Access = private)

        function scaledSequence = scalePilot(obj, targetPower)
            baseSeq = obj.Config.PilotSequence;

            if upper(obj.Config.PilotType) == "SINGLE"
                scaledSequence = sqrt(targetPower) * (baseSeq / abs(baseSeq));
            elseif upper(obj.Config.PilotType) == "CAZAC"
                scaledSequence = sqrt(targetPower) * baseSeq;
            else
                scaledSequence = baseSeq;
            end

        end

        function txSignal = addPrefix(obj, timeFrame)
            prefixLength = obj.Config.PrefixLength;
            N = obj.Config.NumDataSubcarriers;

            if upper(obj.Config.WaveformType) == "AFDM"
                % AFDM: Chirp Periodic Prefix (CPP)
                phaseVector = exp(-1j * 2 * pi * obj.Config.ChirpParam1 * ...
                    (N ^ 2 + 2 * N * (-prefixLength:-1).'));
                preFix = timeFrame(end - prefixLength + 1:end) .* phaseVector;
            else
                % OFDM: 标准 Cyclic Prefix (CP)
                preFix = timeFrame(end - prefixLength + 1:end);
            end

            txSignal = [preFix; timeFrame];
        end

    end

end
