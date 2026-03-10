classdef GiFreeSystem < handle
    % GiFreeSystem  GI-Free AFDM 系统级编排器 (handle 类)
    %
    %   作为顶层 handle 对象, 持有并组合所有子模块:
    %     Transmitter → ChannelBuilder → Receiver (内含 Estimator)
    %
    %   所有子模块共享同一个 GiFreeConfig handle, 修改配置后
    %   无需重新构造任何模块, 新参数自动生效.
    %
    %   用法:
    %     cfg = GiFreeConfig();
    %     cfg.NumSubcarriers = 512;   % 修改参数
    %     sys = GiFreeSystem(cfg);     % 一次构造
    %     result = sys.runTrial(snrLin, noisePow);

    properties (SetAccess = private)
        Config % GiFreeConfig handle
        Transmitter % GiFreeTransmitter handle
        ChannelBuilder % FractionalChannelBuilder handle
        Estimator % GiFreeEstimator handle
        Receiver % GiFreeReceiver handle
    end

    methods

        %% ---------- 构造函数: 一站式组装 ----------
        function obj = GiFreeSystem(cfg)
            % GiFreeSystem  根据配置对象构造完整系统
            %   自动创建并连接所有子模块, 共享同一个 Config handle.
            arguments
                cfg (1, 1) GiFreeConfig
            end

            obj.Config = cfg;
            obj.Transmitter = GiFreeTransmitter(cfg);
            obj.ChannelBuilder = FractionalChannelBuilder(cfg);
            obj.Estimator = GiFreeEstimator(cfg, obj.ChannelBuilder);
            obj.Receiver = GiFreeReceiver(cfg, obj.Estimator);
        end

        %% ========== 单次试验 ==========
        function trialResult = runTrial(obj, dataSnrLinear, noisePower)
            % runTrial  执行一次完整的 GI-Free 系统仿真
            %
            %   输出结构体包含:
            %     bitErrorsSys : 系统比特错误数
            %     bitErrorsRef : Perfect CSI 比特错误数
            %     totalBits    : 总比特数
            %     mseSystem    : 信道估计 NMSE

            numSc = obj.Config.NumSubcarriers;
            bitsPerSymbol = log2(obj.Config.ModulationOrder);

            % ---- 发射 ----
            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLinear);

            % ---- 信道 ----
            [~, hEffTrue] = obj.ChannelBuilder.generateChannel();
            noiseVec = sqrt(noisePower / 2) * ...
                (randn(numSc, 1) + 1j * randn(numSc, 1));
            rxSignal = hEffTrue * txFrame + noiseVec;

            % ---- 接收 ----
            [detIndicesSys, mseSys, ~] = obj.Receiver.receive( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower);

            % ---- Perfect CSI 对比 ----
            detIndicesRef = obj.perfectCsiDetect( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower);

            % ---- 评估 ----
            [bitErrSys, totalBits] = GiFreeSystem.computeBer( ...
                txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~] = GiFreeSystem.computeBer( ...
                txDataIndices, detIndicesRef, bitsPerSymbol);

            trialResult.bitErrorsSys = bitErrSys;
            trialResult.bitErrorsRef = bitErrRef;
            trialResult.totalBits = totalBits;
            trialResult.mseSystem = mseSys;
        end

        %% ========== Perfect CSI 检测 ==========
        function detectedIndices = perfectCsiDetect(obj, ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower)
            % perfectCsiDetect  已知真实信道的 LMMSE 数据检测 (性能上界)

            numSc = obj.Config.NumSubcarriers;
            modOrder = obj.Config.ModulationOrder;

            pilotFrame = zeros(numSc, 1);
            pilotFrame(1) = obj.Config.PilotAmplitude;
            cleanDataSig = rxSignal - hEffTrue * pilotFrame;

            regParam = noisePower / dataSnrLinear;
            hData = hEffTrue(:, 2:numSc);
            estData = (hData' * hData + regParam * speye(numSc - 1)) \ ...
                (hData' * cleanDataSig);

            detectedIndices = qamdemod( ...
                estData / sqrt(dataSnrLinear), modOrder, 'UnitAveragePower', true);
        end

    end

    %% ========== 纯数学工具 — 保留 Static ==========
    methods (Static)

        function [bitErrors, totalBits] = computeBer(txIndices, rxIndices, bitsPerSymbol)
            % computeBer  比较发射与接收符号索引, 统计比特错误
            txBits = de2bi(txIndices, bitsPerSymbol, 'left-msb');
            rxBits = de2bi(rxIndices, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
        end

    end

end
