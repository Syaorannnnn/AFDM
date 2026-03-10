classdef GiFreeTransmitter < handle
    % GiFreeTransmitter  GI-Free AFDM 发射机 (handle 类)
    %
    %   持有 GiFreeConfig 引用, 在 DAFT 域第 0 号子载波放置导频,
    %   其余子载波承载 QAM 数据符号。
    %
    %   用法:
    %     cfg = GiFreeConfig();
    %     tx  = GiFreeTransmitter(cfg);
    %     [txFrame, txDataIndices] = tx.transmit(dataSnrLinear);

    properties (SetAccess = private)
        Config % GiFreeConfig handle 引用
    end

    methods

        function obj = GiFreeTransmitter(cfg)
            % GiFreeTransmitter  构造函数
            arguments
                cfg (1, 1) GiFreeConfig
            end

            obj.Config = cfg;
        end

        function [txFrame, txDataIndices] = transmit(obj, dataSnrLinear)
            % transmit  生成一帧发射信号
            %   txFrame       : N x 1 复数发射向量
            %   txDataIndices : (N-1) x 1 QAM 符号索引 (用于后续 BER 计算)
            numSc    = obj.Config.NumSubcarriers;
            modOrder = obj.Config.ModulationOrder;

              % 随机数据符号
            txDataIndices = randi([0, modOrder - 1], numSc - 1, 1);
            qamSymbols    = qammod(txDataIndices, modOrder, 'UnitAveragePower', true);

              % 组帧: [导频; 数据]
            txFrame        = zeros(numSc, 1);
            txFrame(1)     = obj.Config.PilotAmplitude;
            txFrame(2:end) = qamSymbols * sqrt(dataSnrLinear);
        end

    end

end
