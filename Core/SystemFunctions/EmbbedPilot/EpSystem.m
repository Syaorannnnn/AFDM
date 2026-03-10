classdef EpSystem < handle
    % EpSystem  EP-AFDM 系统级编排器 (单导频版)
    %
    % 用法:
    %   cfg = EpConfig();
    %   cfg.TotalSubcarriers = 256;  % 空口总采样 (与 GiFreeConfig.NumSubcarriers 对齐)
    %   sys = EpSystem(cfg);
    %   result = sys.runTrial(snrDb, pathDelays, pathDopplers, pathGains);

    properties (SetAccess = private)
        Config
        Transmitter
        Receiver
    end

    methods

        function obj = EpSystem(cfg)
            arguments
                cfg (1, 1) EpConfig
            end

            obj.Config = cfg;
            obj.Transmitter = EpTransmitter(cfg);
            obj.Receiver = EpReceiver(cfg);
        end

        %% ========== 单次试验 (估计 CSI) ==========
        function result = runTrial(obj, snrDb, pathDelays, pathDopplers, pathGains)
            totalSc = obj.Config.TotalSubcarriers;
            bps = obj.Config.BitsPerSymbol;

            snrLin = 10 ^ (snrDb / 10);
            pilotPower = 10 ^ (obj.Config.PilotSnr / 10);
            noisePower = 1 / snrLin;

            [txSig, txData] = obj.Transmitter.transmit(pilotPower);

            physH = LtvChannel(totalSc, pathDelays, pathDopplers, pathGains);

            noiseVec = sqrt(noisePower / 2) * (randn(totalSc, 1) + 1j * randn(totalSc, 1));
            rxSig = physH * txSig + noiseVec;

            obj.Receiver.CsiMode = "Estimated";
            rxData = obj.Receiver.receive(rxSig, noisePower, physH, pilotPower);

            [result.bitErrors, result.totalBits] = EpSystem.computeBer(txData, rxData, bps);
        end

        %% ========== 单次试验 (完美 CSI) ==========
        function result = runTrialPerfectCsi(obj, snrDb, pathDelays, pathDopplers, pathGains)
            totalSc = obj.Config.TotalSubcarriers;
            bps = obj.Config.BitsPerSymbol;

            snrLin = 10 ^ (snrDb / 10);
            pilotPower = 10 ^ (obj.Config.PilotSnr / 10);
            noisePower = 1 / snrLin;

            [txSig, txData] = obj.Transmitter.transmit(pilotPower);

            physH = LtvChannel(totalSc, pathDelays, pathDopplers, pathGains);

            noiseVec = sqrt(noisePower / 2) * (randn(totalSc, 1) + 1j * randn(totalSc, 1));
            rxSig = physH * txSig + noiseVec;

            obj.Receiver.CsiMode = "Perfect";
            rxData = obj.Receiver.receive(rxSig, noisePower, physH, pilotPower);

            [result.bitErrors, result.totalBits] = EpSystem.computeBer(txData, rxData, bps);
        end

    end

    methods (Static)

        function [bitErrors, totalBits] = computeBer(txIdx, rxIdx, bitsPerSymbol)
            txBits = de2bi(txIdx, bitsPerSymbol, 'left-msb');
            rxBits = de2bi(rxIdx, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
        end

    end

end