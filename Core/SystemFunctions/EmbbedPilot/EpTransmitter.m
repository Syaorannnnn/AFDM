classdef EpTransmitter < handle
    % EpTransmitter: EP-AFDM 发射机 (单导频版)
    %
    % 信号链: QAM 调制 → 单导频插入 → IDAFT → 添加 CPP

    properties (Access = private)
        Config
    end

    methods (Access = public)

        function obj = EpTransmitter(configObj)
            obj.Config = configObj;
        end

        function [txSignal, originalData] = transmit(obj, pilotPower)
            % transmit  生成一帧发射信号
            %   pilotPower : 导频符号功率 (线性)
            %   txSignal   : TotalSubcarriers × 1 (含 CPP)
            %   originalData : NumActiveCarriers × 1 QAM 符号索引

            N = obj.Config.NumDataSubcarriers;

            % 随机数据 → QAM
            originalData = randi([0, obj.Config.ModulationOrder - 1], ...
                obj.Config.NumActiveCarriers, 1);
            qamSymbols = qammod(originalData, obj.Config.ModulationOrder, ...
                'UnitAveragePower', true);

            % 组装 DAFT 域帧
            daftFrame = zeros(N, 1);

            % 插入单导频 (缩放到指定功率)
            daftFrame(obj.Config.PilotIndex) = sqrt(pilotPower);

            % 映射数据
            daftFrame(obj.Config.ActiveIndices) = qamSymbols;

            % IDAFT 变换到时域
            if upper(obj.Config.WaveformType) == "AFDM"
                timeFrame = AfdmTransforms.idaft(daftFrame, ...
                    obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                timeFrame = AfdmTransforms.idft(daftFrame);
            end

            % 添加 CPP 前缀
            txSignal = obj.addPrefix(timeFrame);
        end

    end

    methods (Access = private)

        function txSignal = addPrefix(obj, timeFrame)
            prefixLen = obj.Config.PrefixLength;
            N = obj.Config.NumDataSubcarriers;

            if upper(obj.Config.WaveformType) == "AFDM"
                phaseVec = exp(-1j * 2 * pi * obj.Config.ChirpParam1 * ...
                    (N ^ 2 + 2 * N * (-prefixLen:-1).'));
                prefix = timeFrame(end - prefixLen + 1:end) .* phaseVec;
            else
                prefix = timeFrame(end - prefixLen + 1:end);
            end

            txSignal = [prefix; timeFrame];
        end

    end

end