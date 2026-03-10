classdef GiFreeReceiver < handle
    % GiFreeReceiver  最终版 GI-Free AFDM Turbo DD-CE 接收机 (handle 类)
    %
    %   持有 GiFreeConfig 与 GiFreeEstimator 的 handle 引用,
    %   实现三阶段接收架构: Phase 1 粗调 → Phase 2 DD Turbo → 后判决精炼
    %
    %   用法:
    %     cfg       = GiFreeConfig();
    %     chBuilder = FractionalChannelBuilder(cfg);
    %     estimator = GiFreeEstimator(cfg, chBuilder);
    %     receiver  = GiFreeReceiver(cfg, estimator);
    %     [detIdx, nmse, hEst] = receiver.receive(rxSignal, hTrue, snrLin, noisePow);
    
    properties (SetAccess = private)
        Config     % GiFreeConfig handle 引用
        Estimator  % GiFreeEstimator handle 引用
    end
    
    methods
        
        %% ---------- 构造函数 ----------
        function obj = GiFreeReceiver(cfg, estimator)
            arguments
                cfg       (1,1) GiFreeConfig
                estimator (1,1) GiFreeEstimator
            end
            obj.Config    = cfg;
            obj.Estimator = estimator;
        end
        
        %% ========== 主接收入口 ==========
        function [detectedIndices, normalizedMse, estEffChannel] = receive( ...
                obj, rxSignal, trueEffChannel, dataSnrLinear, noisePower)
            
            numSc      = obj.Config.NumSubcarriers;
            modOrder   = obj.Config.ModulationOrder;
            totalIter  = obj.Config.MaxSicIterations;
            pilotAmp   = obj.Config.PilotAmplitude;
            
            dataIndices = (2:numSc)';
            numData     = numSc - 1;
            pilotFrame  = zeros(numSc, 1);
            pilotFrame(1) = pilotAmp;
            
            % 预计算星座图
            qamConstellation    = qammod((0:modOrder-1)', modOrder, 'UnitAveragePower', true);
            scaledConstellation = qamConstellation * sqrt(dataSnrLinear);
            sparseEye           = speye(numData);
            
            % 理论最优正则化
            perfectCsiReg = noisePower / dataSnrLinear;
            
            % ---- 先验干扰功率计算 ----
            dataInterfPower     = (numData / numSc) * dataSnrLinear;
            effectiveNoisePower = noisePower + dataInterfPower;
            regParam            = effectiveNoisePower / dataSnrLinear;
            ompRegParam         = dataInterfPower / (pilotAmp^2);
            
            % ---- Phase 分配 ----
            numPhase1       = max(round(totalIter * 0.30), 2);
            numPostDecision = 1;
            numPhase2       = max(totalIter - numPhase1 - numPostDecision, 3);
            
            % ---- 初始 OMP + LMMSE ----
            [estPaths, hEff] = obj.Estimator.estimateByOmp(rxSignal, ompRegParam);
            
            cleanDataSig = rxSignal - hEff * pilotFrame;
            estFullSig   = GiFreeReceiver.lmmseDetect( ...
                cleanDataSig, hEff, regParam, numSc, sparseEye);
            
            % ================================================================
            %   Phase 1: 干扰感知粗调
            % ================================================================
            for iter = 1:numPhase1
                dataFrame   = zeros(numSc, 1);
                rawDataEst  = estFullSig(dataIndices);
                softSymbols = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, effectiveNoisePower);
                dataFrame(dataIndices) = softSymbols;
                
                cleanPilotSig = rxSignal - hEff * dataFrame;
                
                pilotResidual     = cleanPilotSig - hEff * pilotFrame;
                residualPower     = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf    = max(residualPower - noisePower, 0);
                effectiveNoisePower = noisePower + residualInterf;
                regParam    = max(effectiveNoisePower / dataSnrLinear, perfectCsiReg);
                ompRegParam = max(residualInterf / (pilotAmp^2), 0);
                
                [estPaths, hEff] = obj.Estimator.estimateByOmp(cleanPilotSig, ompRegParam);
                
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end
            
            % ================================================================
            %   Phase 2: DD Turbo 精调
            % ================================================================
            for iter = 1:numPhase2
                
                rawDataEst     = estFullSig(dataIndices);
                outputNoiseVar = GiFreeReceiver.estimateOutputNoise( ...
                    rawDataEst, scaledConstellation, noisePower);
                softSymbols    = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar);
                
                txFrameEst    = zeros(numSc, 1);
                txFrameEst(1) = pilotAmp;
                txFrameEst(dataIndices) = softSymbols;
                
                % ID2P 消除 → OMP 路径重检测
                dataFrameDd = zeros(numSc, 1);
                dataFrameDd(dataIndices) = softSymbols;
                cleanPilotSig = rxSignal - hEff * dataFrameDd;
                
                pilotResidual  = cleanPilotSig - hEff * pilotFrame;
                residualPower  = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf = max(residualPower - noisePower, 0);
                ompRegParam    = max(residualInterf / (pilotAmp^2), 0);
                
                [estPaths, ~] = obj.Estimator.estimateByOmp(cleanPilotSig, ompRegParam);
                
                % DD 精炼
                txFramePower = real(txFrameEst' * txFrameEst) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);
                
                if iter <= 3
                    [estPaths, hEff] = obj.Estimator.refineDopplerByDd( ...
                        estPaths, rxSignal, txFrameEst, ddRegParam);
                else
                    [estPaths, hEff] = obj.Estimator.refineGainByDd( ...
                        estPaths, rxSignal, txFrameEst, ddRegParam);
                end
                
                % 截断均值正则化
                trimAlpha         = max(0.15 - 0.02 * iter, 0.05);
                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameEst, trimAlpha);
                
                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);
                % effectiveNoisePower = max(trimmedNoisePower, noisePower);
                
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end
            
            % ================================================================
            %   后判决精炼
            % ================================================================
            for iter = 1:numPostDecision
                normalizedEst = estFullSig(dataIndices) / sqrt(dataSnrLinear);
                hardDecIdx    = qamdemod(normalizedEst, modOrder, 'UnitAveragePower', true);
                hardDecSym    = qammod(hardDecIdx, modOrder, 'UnitAveragePower', true) * sqrt(dataSnrLinear);
                
                txFrameHard    = zeros(numSc, 1);
                txFrameHard(1) = pilotAmp;
                txFrameHard(dataIndices) = hardDecSym;
                
                txFramePower = real(txFrameHard' * txFrameHard) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);
                
                [estPaths, hEff] = obj.Estimator.refineGainByDd( ...
                    estPaths, rxSignal, txFrameHard, ddRegParam);
                
                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameHard, 0.03);
                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);
                
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end
            
            % ---- 最终解调 ----
            normalizedDataSym = estFullSig(dataIndices) / sqrt(dataSnrLinear);
            detectedIndices   = qamdemod(normalizedDataSym, modOrder, 'UnitAveragePower', true);
            estEffChannel     = hEff;
            
            normalizedMse = norm(full(hEff) - full(trueEffChannel), 'fro')^2 / ...
                max(norm(full(trueEffChannel), 'fro')^2, 1e-20);
        end
        
    end
    
    %% ========== 无状态数学工具 — 保留 Static ==========
    methods (Static)
        
        function estSig = lmmseDetect(cleanSig, hEff, regParam, numSc, sparseEye)
            % lmmseDetect  LMMSE 均衡 (纯数学, 不依赖任何对象状态)
            hData   = hEff(:, 2:numSc);
            estData = (hData' * hData + regParam * sparseEye) \ (hData' * cleanSig);
            estSig          = zeros(numSc, 1);
            estSig(2:numSc) = estData;
        end
        
        function softSym = computeSoftSymbols(rawEst, constellation, noiseVar)
            % computeSoftSymbols  MMSE 软符号 (纯数学)
            numSymbols = length(rawEst);
            softSym    = zeros(numSymbols, 1);
            sigma2     = noiseVar;
            
            for k = 1:numSymbols
                distances  = abs(rawEst(k) - constellation).^2;
                logWeights = -distances / sigma2;
                logWeights = logWeights - max(logWeights);
                weights    = exp(logWeights);
                weights    = weights / sum(weights);
                softSym(k) = sum(weights .* constellation);
            end
        end
        
        function outputNoiseVar = estimateOutputNoise(rawEst, scaledConstellation, noisePower)
            % estimateOutputNoise  输出 SINR 噪声估计 (纯数学)
            numSymbols    = length(rawEst);
            distToNearest = zeros(numSymbols, 1);
            
            for k = 1:numSymbols
                allDist = abs(rawEst(k) - scaledConstellation).^2;
                distToNearest(k) = min(allDist);
            end
            
            sortedDist     = sort(distToNearest, 'ascend');
            numKeep        = max(floor(numSymbols * 0.9), 1);
            outputNoiseVar = max(mean(sortedDist(1:numKeep)), noisePower);
        end
        
    end
end
