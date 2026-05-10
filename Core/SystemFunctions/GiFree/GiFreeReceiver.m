classdef GiFreeReceiver < handle
% GIFREERECEIVER - CLIP 三阶段接收处理主流程。
    properties (SetAccess = private)
        Config
        Estimator
    end

    methods

        % GIFREERECEIVER 构造函数，绑定配置与估计器实例。
        function obj = GiFreeReceiver(cfg, estimator)
            arguments
                cfg       (1,1) GiFreeConfig
                estimator (1,1) GiFreeEstimator
            end
            obj.Config    = cfg;
            obj.Estimator = estimator;
        end

        % RECEIVE 执行 CLIP 接收链路并输出检测结果与信道估计误差。
        function [detectedIndices, normalizedMse, estEffChannel, rxDiag] = receive( ...
                obj, rxSignal, trueEffChannel, dataSnrLin, noisePowerLin)

            numSc      = obj.Config.NumSubcarriers;
            modOrder   = obj.Config.ModulationOrder;
            totalIter  = obj.Config.MaxSicIterations;
            numData    = obj.Config.NumDataSymbols;

            % NOTE: 配置对象返回的 PilotPos1/DataPos1 已是 MATLAB 1-based 索引。
            pilotPos1  = obj.Config.PilotPos1;
            dataPos1   = obj.Config.DataPos1;
            pilotSeq   = obj.Config.PilotSequence;
            ampPP      = obj.Config.PerPilotAmplitude;
            pilotAmpTot = obj.Config.PilotAmplitude;

            pilotFrame = zeros(numSc, 1);
            pilotFrame(pilotPos1) = ampPP * pilotSeq;

            qamConstellation    = qammod((0:modOrder-1)', modOrder, 'UnitAveragePower', true);
            scaledConstellation = qamConstellation * sqrt(dataSnrLin);

            % 正则化参数
            perfectCsiReg       = noisePowerLin / dataSnrLin;
            dataInterfPower     = (numData / numSc) * dataSnrLin;
            effectiveNoisePower = noisePowerLin + dataInterfPower;
            regParam            = effectiveNoisePower / dataSnrLin;
            ompRegParam         = min(dataInterfPower / (pilotAmpTot^2), obj.Config.OmpRegUpperBound);
            baseSearchMode      = lower(obj.Config.EarlyPathSearchMode);
            rxDiag              = GiFreeReceiver.createEmptyRxDiag();
            prevOmpDiag         = [];
            lastOmpDiag         = [];
            lastOmpSupportChanged = strcmpi(baseSearchMode, 'beam');

            prioritizePhase1Budget = obj.Config.EnableProgressiveCfar && ...
                obj.Config.EnablePathStabilityGate;
            [numPhase1, numPhase2, numPostDecision] = ...
                GiFreeReceiver.splitPhaseBudget( ...
                    totalIter, obj.Config.EnableDfp, prioritizePhase1Budget, ...
                    obj.Config.Phase1BudgetRatio, obj.Config.Phase1BudgetRatioAlt);
            rxDiag.phase1IterationBudget = numPhase1;
            rxDiag.phase2IterationBudget = numPhase2;
            rxDiag.postDecisionIterationBudget = numPostDecision;

            if obj.Config.EnablePhase0Bootstrap && obj.Config.DirichletRadius > 0
                [phase1InputSignal, ~, bootstrapDiag] = obj.bootstrapPhase0( ...
                    rxSignal, ompRegParam, regParam, ...
                    numSc, dataPos1, scaledConstellation, effectiveNoisePower, pilotFrame);
                rxDiag.bootstrapDiag = bootstrapDiag;
            else
                phase1InputSignal = rxSignal;
            end

            [estPaths, effectiveChannel, ompDiag] = obj.Estimator.estimateByOmp( ...
                phase1InputSignal, ompRegParam, ones(numSc, 1), baseSearchMode);
            rxDiag = GiFreeReceiver.updateOmpCounters(rxDiag, ompDiag);
            rxDiag.phase1ResidualNorm(end+1, 1) = GiFreeReceiver.extractFinalResidualNorm(ompDiag);
            lastOmpDiag = ompDiag;
            pathHistory = GiFreeReceiver.initializePathHistory(estPaths);

            cleanDataSig = rxSignal - effectiveChannel * pilotFrame;
            estFullSig   = GiFreeReceiver.lmmseDetect( ...
                cleanDataSig, effectiveChannel, regParam, numSc, dataPos1);

            for iter = 1:numPhase1
                stablePaths = GiFreeReceiver.extractStablePaths( ...
                    pathHistory, obj.Config.PathStabilityThreshold, ...
                    obj.Config.EnablePathStabilityGate);
                if isempty(stablePaths)
                    cleaningChannel = effectiveChannel;
                else
                    cleaningChannel = obj.Estimator.ChannelBuilder.buildEffectiveChannel(stablePaths);
                end

                rawDataEst  = estFullSig(dataPos1);
                [~, confidence] = GiFreeReceiver.computeFeedbackSymbols( ...
                    rawDataEst, scaledConstellation, effectiveNoisePower, obj.Config.FeedbackMode);
                cleaningSymbols = GiFreeReceiver.computeCleaningSymbols( ...
                    rawDataEst, scaledConstellation, effectiveNoisePower, obj.Config.CleaningFeedbackMode);
                keepRatio = GiFreeReceiver.phaseKeepRatio( ...
                    iter, numPhase1, 'phase1', obj.Config.EnableConfidenceGating, ...
                    obj.Config.CleaningGateMode, confidence, ...
                    obj.Config.AdaptConfScaleDiv, obj.Config.AdaptConfScaleMin, obj.Config.AdaptConfScaleMax, obj.Config.AdaptConfKeepCap);
                gatedCleaning = GiFreeReceiver.applyConfidenceGate(cleaningSymbols, confidence, keepRatio);

                cleanPilotSig = rxSignal - cleaningChannel * ...
                    GiFreeReceiver.buildDataFrame(numSc, dataPos1, gatedCleaning);
                obsWeights = GiFreeReceiver.buildObservationWeights( ...
                    cleaningChannel, numSc, dataPos1, cleaningSymbols, gatedCleaning, ...
                    obj.Config.UseWeightedPilotMetric);

                pilotResidual       = cleanPilotSig - cleaningChannel * pilotFrame;
                residualPower       = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf      = max(residualPower - noisePowerLin, 0);
                effectiveNoisePower = noisePowerLin + residualInterf;
                regParam    = max(effectiveNoisePower / dataSnrLin, perfectCsiReg);
                phase1BaseOmpRegParam = max(residualInterf / (pilotAmpTot^2), 0);
                effectiveDetectionStrictness = GiFreeReceiver.computePhase1DetectionStrictness( ...
                    iter, numPhase1, obj.Config.EnableProgressiveCfar, ...
                    obj.Config.ProgressiveCfarInitScale, obj.Config.ProgressiveCfarFinalScale);
                ompRegParam = min(phase1BaseOmpRegParam * effectiveDetectionStrictness, obj.Config.Phase1OmpRegUpperBound);

                prevPaths = estPaths;
                searchMode = GiFreeReceiver.selectPhase1OmpMode( ...
                    baseSearchMode, obj.Config.BeamStrategyMode, obj.Config.DirichletRadius, ...
                    obj.Config.BeamAdaptiveSpreadThreshold, obj.Config.BeamAdaptiveImproveRatio, ...
                    iter, prevOmpDiag, lastOmpDiag, lastOmpSupportChanged);
                [estPaths, effectiveChannel, ompDiag] = obj.Estimator.estimateByOmp( ...
                    cleanPilotSig, ompRegParam, obsWeights, searchMode);
                supportChanged = GiFreeReceiver.hasSupportChange(prevPaths, estPaths);
                if supportChanged
                    rxDiag.phase1SupportChangeCount = rxDiag.phase1SupportChangeCount + 1;
                end
                rxDiag = GiFreeReceiver.updateOmpCounters(rxDiag, ompDiag);
                rxDiag.phase1ResidualNorm(end+1, 1) = GiFreeReceiver.extractFinalResidualNorm(ompDiag);
                prevOmpDiag = lastOmpDiag;
                lastOmpDiag = ompDiag;
                lastOmpSupportChanged = supportChanged;
                pathHistory = GiFreeReceiver.updatePathHistory( ...
                    pathHistory, estPaths, obj.Config.PathStabilityDopplerTolerance);
                stablePaths = GiFreeReceiver.extractStablePaths( ...
                    pathHistory, obj.Config.PathStabilityThreshold, ...
                    obj.Config.EnablePathStabilityGate);
                rxDiag.phase1StablePathCount(end+1, 1) = size(stablePaths, 1);
                rxDiag.phase1TotalPathCount(end+1, 1) = size(estPaths, 1);
                rxDiag.phase1EffectiveDetectionStrictness(end+1, 1) = effectiveDetectionStrictness;

                cleanDataSig = rxSignal - effectiveChannel * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, effectiveChannel, regParam, numSc, dataPos1);
            end

            stagnationCount = 0;
            % 用 Phase1 末尾的 effectiveNoisePower 初始化软判决置信度下界。
            % 后续每轮迭代由 trimmedNoisePower 更新：H 越差残差越大，下界越高，
            % 软判决越保守，从而形成自校正的负反馈机制，抑制高 SNR 下的 Turbo 发散。
            softDecNoiseFloor = effectiveNoisePower;
            % NOTE: 迭代精炼阶段交替执行支撑更新与数据导向精炼。
            for iter = 1:numPhase2
                rawDataEst     = estFullSig(dataPos1);
                outputNoiseVar = GiFreeReceiver.selectOutputNoiseVar( ...
                    rawDataEst, scaledConstellation, softDecNoiseFloor, obj.Config.UseAdaptiveOutputNoise, ...
                    obj.Config.OutputNoiseKeepRatio);
                [softSymbols, confidence] = GiFreeReceiver.computeFeedbackSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar, obj.Config.FeedbackMode);
                cleaningSymbols = GiFreeReceiver.computeCleaningSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar, obj.Config.CleaningFeedbackMode);
                keepRatio = GiFreeReceiver.phaseKeepRatio( ...
                    iter, numPhase2, 'phase2', obj.Config.EnableConfidenceGating, ...
                    obj.Config.CleaningGateMode, confidence, ...
                    obj.Config.AdaptConfScaleDiv, obj.Config.AdaptConfScaleMin, obj.Config.AdaptConfScaleMax, obj.Config.AdaptConfKeepCap);
                gatedCleaning = GiFreeReceiver.applyConfidenceGate(cleaningSymbols, confidence, keepRatio);

                txFrameEst = GiFreeReceiver.buildTxFrame( ...
                    numSc, pilotPos1, dataPos1, ampPP, pilotSeq, softSymbols);
                prevResNorm = norm(rxSignal - effectiveChannel * txFrameEst)^2 / numSc;

                cleanPilotSig = rxSignal - effectiveChannel * ...
                    GiFreeReceiver.buildDataFrame(numSc, dataPos1, gatedCleaning);
                obsWeights = GiFreeReceiver.buildObservationWeights( ...
                    effectiveChannel, numSc, dataPos1, cleaningSymbols, gatedCleaning, ...
                    obj.Config.UseWeightedPilotMetric);

                pilotResidual  = cleanPilotSig - effectiveChannel * pilotFrame;
                residualPower  = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf = max(residualPower - noisePowerLin, 0);
                ompRegParam    = min(max(residualInterf / (pilotAmpTot^2), 0), obj.Config.OmpRegUpperBound);

                prevPaths = estPaths;
                searchMode = GiFreeReceiver.selectPhase2OmpMode( ...
                    baseSearchMode, lastOmpSupportChanged, stagnationCount, lastOmpDiag);
                [ompEstPaths, hEffOmp, ompDiag] = obj.Estimator.estimateByOmp( ...
                    cleanPilotSig, ompRegParam, obsWeights, searchMode);
                supportChanged = GiFreeReceiver.hasSupportChange(prevPaths, ompEstPaths);
                rxDiag = GiFreeReceiver.updateOmpCounters(rxDiag, ompDiag);
                prevOmpDiag = lastOmpDiag;
                lastOmpDiag = ompDiag;

                txFramePower = real(txFrameEst' * txFrameEst) / numSc;
                ddRegParam   = max(noisePowerLin / txFramePower, 1e-6);
                useDopplerRefine = obj.Config.DirichletRadius > 0 && ...
                    (iter == 1 || supportChanged || stagnationCount >= 1 || mod(iter, 2) == 1);

                if useDopplerRefine
                    [candEstPaths, candHEff] = obj.Estimator.refineDopplerByDd( ...
                        ompEstPaths, rxSignal, txFrameEst, ddRegParam);
                else
                    [candEstPaths, candHEff] = obj.Estimator.refineGainByDd( ...
                        ompEstPaths, rxSignal, txFrameEst, ddRegParam);
                end

                ompResNorm  = norm(rxSignal - hEffOmp * txFrameEst)^2 / numSc;
                candResNorm = norm(rxSignal - candHEff * txFrameEst)^2 / numSc;
                if candResNorm <= ompResNorm
                    bestNewPaths   = candEstPaths;
                    bestNewHEff    = candHEff;
                    bestNewResNorm = candResNorm;
                else
                    bestNewPaths   = ompEstPaths;
                    bestNewHEff    = hEffOmp;
                    bestNewResNorm = ompResNorm;
                end

                localCleanPilotSig = rxSignal - bestNewHEff * ...
                    GiFreeReceiver.buildDataFrame(numSc, dataPos1, gatedCleaning);
                [injCandidate, ~, ~] = obj.Estimator.proposeResidualPath( ...
                    localCleanPilotSig, bestNewPaths, ompRegParam, obj.Config.ResidualPathThresholdScale, obsWeights);
                if obj.Config.EnableResidualInjection && ~isempty(injCandidate)
                    unionPaths = GiFreeReceiver.pruneWeakPaths( ...
                        [bestNewPaths; injCandidate], obj.Config.NumPathsUpper);
                    if useDopplerRefine
                        [unionPaths, unionHEff] = obj.Estimator.refineDopplerByDd( ...
                            unionPaths, rxSignal, txFrameEst, ddRegParam);
                    else
                        [unionPaths, unionHEff] = obj.Estimator.refineGainByDd( ...
                            unionPaths, rxSignal, txFrameEst, ddRegParam);
                    end
                    unionResNorm = norm(rxSignal - unionHEff * txFrameEst)^2 / numSc;
                    if unionResNorm < bestNewResNorm * obj.Config.UnionImproveRatio
                        bestNewPaths   = unionPaths;
                        bestNewHEff    = unionHEff;
                        bestNewResNorm = unionResNorm;
                        supportChanged = true;
                    end
                end

                acceptTol = obj.Config.Phase2AcceptTol;
                if supportChanged
                    acceptTol = obj.Config.Phase2AcceptTolWithChange;
                end
                improveTol = obj.Config.Phase2ImproveTol;
                acceptedSupportChanged = false;

                if bestNewResNorm <= prevResNorm * (1 + acceptTol)
                    estPaths = bestNewPaths;
                    effectiveChannel = bestNewHEff;
                    acceptedSupportChanged = GiFreeReceiver.hasSupportChange(prevPaths, estPaths);
                    if bestNewResNorm <= prevResNorm * (1 - improveTol)
                        stagnationCount = 0;
                    else
                        stagnationCount = min(stagnationCount + 1, 3);
                    end
                else
                    stagnationCount = min(stagnationCount + 1, 3);
                end

                trimAlpha         = max(obj.Config.TrimAlphaInit - obj.Config.TrimAlphaDecay * iter, obj.Config.TrimAlphaMin);
                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, effectiveChannel, txFrameEst, trimAlpha);
                regParam          = max(trimmedNoisePower / dataSnrLin, perfectCsiReg);
                % 更新软判决置信度下界：本轮全帧残差反映当前信道估计质量。
                % 下轮 selectOutputNoiseVar 将以此为最低参考，防止对错误决策过度置信。
                softDecNoiseFloor = trimmedNoisePower;
                if acceptedSupportChanged
                    rxDiag.phase2SupportChangeCount = rxDiag.phase2SupportChangeCount + 1;
                end
                rxDiag.phase2ResidualNorm(end+1, 1) = ...
                    norm(rxSignal - effectiveChannel * txFrameEst)^2 / numSc;
                lastOmpSupportChanged = supportChanged || acceptedSupportChanged;

                cleanDataSig = rxSignal - effectiveChannel * pilotFrame;
                estFullSig   = SparseLmmseSolver.solve( ...
                    cleanDataSig, effectiveChannel, regParam, numSc, ...
                    'MaxIter', 20, 'Tolerance', 1e-6, ...
                    'WarmStart', estFullSig, 'DataPos1', dataPos1);
            end

            for iter = 1:numPostDecision
                rawDataEst = estFullSig(dataPos1);
                outputNoiseVar = GiFreeReceiver.selectOutputNoiseVar( ...
                    rawDataEst, scaledConstellation, noisePowerLin, obj.Config.UseAdaptiveOutputNoise, ...
                    obj.Config.OutputNoiseKeepRatio);
                [softSymDfp, confidenceDfp] = GiFreeReceiver.computeFeedbackSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar, obj.Config.FeedbackMode);
                cleaningSymbolsDfp = GiFreeReceiver.computeCleaningSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar, obj.Config.CleaningFeedbackMode);
                keepRatio = GiFreeReceiver.phaseKeepRatio( ...
                    iter, numPostDecision, 'dfp', obj.Config.EnableConfidenceGating, ...
                    obj.Config.CleaningGateMode, confidenceDfp, ...
                    obj.Config.AdaptConfScaleDiv, obj.Config.AdaptConfScaleMin, obj.Config.AdaptConfScaleMax, obj.Config.AdaptConfKeepCap);
                gatedCleaningDfp = GiFreeReceiver.applyConfidenceGate( ...
                    cleaningSymbolsDfp, confidenceDfp, keepRatio);

                txFrameSoft = GiFreeReceiver.buildTxFrame( ...
                    numSc, pilotPos1, dataPos1, ampPP, pilotSeq, softSymDfp);
                txFramePower = real(txFrameSoft' * txFrameSoft) / numSc;
                ddRegParam   = max(noisePowerLin / txFramePower, 1e-6);
                [estPaths, effectiveChannel] = obj.Estimator.refineGainByDd( ...
                    estPaths, rxSignal, txFrameSoft, ddRegParam);

                txFrameSoftGated = GiFreeReceiver.buildTxFrame( ...
                    numSc, pilotPos1, dataPos1, ampPP, pilotSeq, gatedCleaningDfp);
                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, effectiveChannel, txFrameSoftGated, 0.03);
                regParam = max(trimmedNoisePower / dataSnrLin, perfectCsiReg);

                cleanDataSig = rxSignal - effectiveChannel * pilotFrame;
                estFullSig   = SparseLmmseSolver.solve( ...
                    cleanDataSig, effectiveChannel, regParam, numSc, ...
                    'MaxIter', 15, 'Tolerance', 1e-6, ...
                    'WarmStart', estFullSig, 'DataPos1', dataPos1);
            end

            % 最终解调
            normalizedDataSym = estFullSig(dataPos1) / sqrt(dataSnrLin);
            detectedIndices   = qamdemod(normalizedDataSym, modOrder, 'UnitAveragePower', true);
            estEffChannel     = effectiveChannel;
            normalizedMse     = norm(full(effectiveChannel) - full(trueEffChannel), 'fro')^2 / ...
                max(norm(full(trueEffChannel), 'fro')^2, 1e-20);

            % 填充数据检测诊断字段
            rxDiag.regParamFinal = regParam;
            txFrameEst = GiFreeReceiver.buildTxFrame( ...
                numSc, pilotPos1, dataPos1, ampPP, pilotSeq, estFullSig(dataPos1));
            dataResidual = rxSignal - effectiveChannel * txFrameEst;
            rxDiag.dataResidualPower = real(dataResidual' * dataResidual) / numSc;
            rxDiag.pilotDataPowerRatio = pilotAmpTot^2 / dataSnrLin;
            rxDiag.finalEstimatedPaths = estPaths;
            rxDiag.finalPathCount = size(estPaths, 1);
            [~, finalConfidence] = GiFreeReceiver.computeFeedbackSymbols( ...
                estFullSig(dataPos1), scaledConstellation, noisePowerLin, obj.Config.FeedbackMode);
            rxDiag.symbolConfidenceMean = mean(finalConfidence);
            rxDiag.symbolConfidenceStd = std(finalConfidence);
            % 导频消除残差诊断
            pilotResidual = rxSignal - effectiveChannel * pilotFrame;
            rxDiag.pilotResidualPower = real(pilotResidual' * pilotResidual) / numSc;
            rxDiag.pilotCancellationErrorNorm = rxDiag.pilotResidualPower / pilotAmpTot^2;
        end

    end

    methods (Access = private)

        function [cleanedSignal, intPaths, bootstrapDiag] = bootstrapPhase0( ...
                obj, rxSignal, ompRegParam, regParam, ...
                numSc, dataPos1, scaledConstellation, effectiveNoisePower, pilotFrame)
            originalDirichletRadius = obj.Config.DirichletRadius;
            cleanupObj = onCleanup(@() GiFreeReceiver.restoreDirichletRadius( ...
                obj.Config, originalDirichletRadius));
            obj.Config.DirichletRadius = 0;

            [intPathsRaw, ~, intOmpDiag] = obj.Estimator.estimateByOmp( ...
                rxSignal, ompRegParam, ones(numSc, 1), 'greedy');

            [intPaths, improveRatioVec] = GiFreeReceiver.pruneBootstrapPathsByResidual( ...
                intPathsRaw, intOmpDiag.weightedResidualNormPerDepth, ...
                obj.Config.MaxBootstrapPaths, obj.Config.BootstrapImproveRatioThreshold);

            obj.Config.DirichletRadius = originalDirichletRadius;

            if isempty(intPaths)
                intChannel = sparse(numSc, numSc);
                localRefineAccepted = false;
                localRefineResGain = 1.0;
            else
                intChannel = obj.Estimator.ChannelBuilder.buildEffectiveChannel(intPaths);
                [intPaths, intChannel, localRefineAccepted, localRefineResGain] = ...
                    obj.refineBootstrapPathsLocally( ...
                        intPaths, intChannel, rxSignal, regParam, effectiveNoisePower, ...
                        numSc, dataPos1, scaledConstellation, pilotFrame);
            end

            cleanDataSig = rxSignal - intChannel * pilotFrame;
            intDataEst = GiFreeReceiver.lmmseDetect( ...
                cleanDataSig, intChannel, regParam, numSc, dataPos1);

            rawDataEst = intDataEst(dataPos1);
            [cleaningSymbols, confidence] = GiFreeReceiver.computeFeedbackSymbols( ...
                rawDataEst, scaledConstellation, effectiveNoisePower, 'soft');

            cleanRatio = min(max(obj.Config.Phase0CleanRatio, 0), 1);
            gatedCleaning = GiFreeReceiver.applyConfidenceGate( ...
                cleaningSymbols, confidence, cleanRatio);

            dataFrame = GiFreeReceiver.buildDataFrame(numSc, dataPos1, gatedCleaning);
            cleanedSignal = rxSignal - intChannel * dataFrame;

            bootstrapDiag = struct();
            bootstrapDiag.numIntPathsDetectedRaw = size(intPathsRaw, 1);
            bootstrapDiag.numIntPathsDetected = size(intPaths, 1);
            bootstrapDiag.intOmpDiag = intOmpDiag;
            bootstrapDiag.bootstrapImproveRatioVec = improveRatioVec;
            bootstrapDiag.localRefineAccepted = localRefineAccepted;
            bootstrapDiag.localRefineResGain = localRefineResGain;
            if isempty(confidence)
                bootstrapDiag.meanConfidence = 0;
            else
                bootstrapDiag.meanConfidence = mean(confidence);
            end
            bootstrapDiag.numSymbolsCleaned = sum(abs(gatedCleaning) > 0);
            clear cleanupObj;
        end

        function [bestPaths, bestChannel, localRefineAccepted, localRefineResGain] = ...
                refineBootstrapPathsLocally( ...
                obj, intPaths, intChannel, rxSignal, regParam, effectiveNoisePower, ...
                numSc, dataPos1, scaledConstellation, pilotFrame)
            bestPaths = intPaths;
            bestChannel = intChannel;
            localRefineAccepted = false;
            localRefineResGain = 1.0;

            if ~obj.Config.EnablePhase0LocalRefine || isempty(intPaths) || obj.Config.DirichletRadius == 0
                return;
            end

            cleanDataSig = rxSignal - intChannel * pilotFrame;
            intDataEst = GiFreeReceiver.lmmseDetect( ...
                cleanDataSig, intChannel, regParam, numSc, dataPos1);
            rawDataEst = intDataEst(dataPos1);
            [bootstrapSoftSymbols, ~] = GiFreeReceiver.computeFeedbackSymbols( ...
                rawDataEst, scaledConstellation, effectiveNoisePower, 'soft');

            txFrameBootstrap = pilotFrame;
            txFrameBootstrap(dataPos1) = bootstrapSoftSymbols;
            txFramePower = real(txFrameBootstrap' * txFrameBootstrap) / numSc;
            ddRegParam = max(effectiveNoisePower / max(txFramePower, 1e-12), 1e-6);

            [refinedPaths, refinedChannel] = obj.Estimator.refineDopplerByDd( ...
                intPaths, rxSignal, txFrameBootstrap, ddRegParam);

            intResNorm = norm(rxSignal - intChannel * txFrameBootstrap)^2 / numSc;
            refinedResNorm = norm(rxSignal - refinedChannel * txFrameBootstrap)^2 / numSc;
            localRefineResGain = intResNorm / max(refinedResNorm, 1e-20);

            if refinedResNorm < intResNorm
                bestPaths = refinedPaths;
                bestChannel = refinedChannel;
                localRefineAccepted = true;
            end
        end

    end

    methods (Static)

        % CREATEEMPTYRXDIAG 构造字段完整的接收机诊断结构体（含数据检测诊断）。
        function rxDiag = createEmptyRxDiag()
            rxDiag = struct( ...
                'phase1SupportChangeCount', 0, ...
                'phase2SupportChangeCount', 0, ...
                'phase1ResidualNorm', zeros(0, 1), ...
                'phase2ResidualNorm', zeros(0, 1), ...
                'phase1StablePathCount', zeros(0, 1), ...
                'phase1TotalPathCount', zeros(0, 1), ...
                'phase1EffectiveDetectionStrictness', zeros(0, 1), ...
                'phase1IterationBudget', 0, ...
                'phase2IterationBudget', 0, ...
                'postDecisionIterationBudget', 0, ...
                'ompCallCount', 0, ...
                'beamCallCount', 0, ...
                'bootstrapDiag', [], ...
                'finalEstimatedPaths', zeros(0, 3), ...
                'finalPathCount', 0, ...
                'regParamFinal', 0, ...
                'dataResidualPower', 0, ...
                'pilotDataPowerRatio', 0, ...
                'symbolConfidenceMean', 0, ...
                'symbolConfidenceStd', 0, ...
                'pilotResidualPower', 0, ...
                'pilotCancellationErrorNorm', 0);
        end

        % UPDATEOMPCOUNTERS 统计接收机内部 OMP 总调用次数与 Beam 次数。
        function rxDiag = updateOmpCounters(rxDiag, ompDiag)
            rxDiag.ompCallCount = rxDiag.ompCallCount + 1;
            if strcmpi(ompDiag.modeUsed, 'beam')
                rxDiag.beamCallCount = rxDiag.beamCallCount + 1;
            end
        end

        % EXTRACTFINALRESIDUALNORM 提取一次 OMP 调用结束时的最终残差。
        function residualNorm = extractFinalResidualNorm(ompDiag)
            if isempty(ompDiag.weightedResidualNormPerDepth)
                residualNorm = inf;
            else
                residualNorm = ompDiag.weightedResidualNormPerDepth(end);
            end
        end

        % COMPUTEOMPIMPROVERATIO 计算连续两次 OMP 之间的残差改进比例。
        function improveRatio = computeOmpImproveRatio(prevOmpDiag, curOmpDiag)
            if isempty(prevOmpDiag) || isempty(curOmpDiag)
                improveRatio = 0;
                return;
            end

            prevResidual = GiFreeReceiver.extractFinalResidualNorm(prevOmpDiag);
            curResidual = GiFreeReceiver.extractFinalResidualNorm(curOmpDiag);
            improveRatio = max(prevResidual - curResidual, 0) / max(prevResidual, 1e-12);
        end

        % ISSUPPORTSTABLE 判断连续两次 OMP 是否已进入支撑稳定状态。
        function isStable = isSupportStable(prevOmpDiag, curOmpDiag, minImproveRatio)
            if isempty(prevOmpDiag) || isempty(curOmpDiag)
                isStable = false;
                return;
            end

            sameSignature = strcmp(prevOmpDiag.finalSupportSignature, curOmpDiag.finalSupportSignature);
            improveRatio = GiFreeReceiver.computeOmpImproveRatio(prevOmpDiag, curOmpDiag);
            isStable = sameSignature && improveRatio < minImproveRatio;
        end

        % SELECTPHASE1OMPMODE 根据 Phase1 迭代进度选择 Beam 或 Greedy。
        function searchMode = selectPhase1OmpMode( ...
                baseSearchMode, beamStrategyMode, spreadWidth, adaptiveSpreadThreshold, ...
                adaptiveImproveRatio, iter, prevOmpDiag, lastOmpDiag, lastSupportChanged)
            if ~strcmpi(baseSearchMode, 'beam')
                searchMode = baseSearchMode;
                return;
            end
            beamStrategyMode = lower(beamStrategyMode);
            if iter == 1
                searchMode = 'beam';
                return;
            end

            improveRatio = GiFreeReceiver.computeOmpImproveRatio(prevOmpDiag, lastOmpDiag);
            if strcmp(beamStrategyMode, 'performance')
                if iter <= 2
                    searchMode = 'beam';
                else
                    searchMode = 'greedy';
                end
                return;
            end

            if strcmp(beamStrategyMode, 'fastadaptive')
                if iter == 2 && (spreadWidth >= adaptiveSpreadThreshold || ...
                        lastSupportChanged || improveRatio >= adaptiveImproveRatio)
                    searchMode = 'beam';
                else
                    searchMode = 'greedy';
                end
                return;
            end

            if GiFreeReceiver.isSupportStable(prevOmpDiag, lastOmpDiag, 0.01)
                searchMode = 'greedy';
            elseif lastSupportChanged || improveRatio >= adaptiveImproveRatio
                searchMode = 'beam';
            else
                searchMode = 'greedy';
            end
        end

        % SELECTPHASE2OMPMODE 根据上轮支撑变化与停滞状态选择 Beam 或 Greedy。
        function searchMode = selectPhase2OmpMode( ...
                baseSearchMode, supportChanged, stagnationCount, lastOmpDiag)
            if ~strcmpi(baseSearchMode, 'beam')
                searchMode = baseSearchMode;
                return;
            end

            lastUsedBeam = ~isempty(lastOmpDiag) && strcmpi(lastOmpDiag.modeUsed, 'beam');
            if supportChanged
                if lastUsedBeam && stagnationCount == 0
                    searchMode = 'greedy';
                else
                    searchMode = 'beam';
                end
            elseif stagnationCount >= 2
                searchMode = 'beam';
            else
                searchMode = 'greedy';
            end
        end

        % LMMSEDETECT 对数据子载波执行闭式 LMMSE 检测。
        function estSig = lmmseDetect(cleanSig, effectiveChannel, regParam, numSc, dataPos1)
            hData   = effectiveChannel(:, dataPos1);
            numData = length(dataPos1);
            estData = (hData' * hData + regParam * speye(numData)) \ (hData' * cleanSig);
            estSig  = zeros(numSc, 1);
            estSig(dataPos1) = estData;
        end

        % COMPUTESOFTSYMBOLS 基于后验权重计算软符号与置信度。
        function [softSym, confidence] = computeSoftSymbols(rawEst, constellation, noiseVar)
            distances  = abs(rawEst - constellation.').^2;
            logWeights = -distances / max(noiseVar, 1e-12);
            logWeights = logWeights - max(logWeights, [], 2);
            weights    = exp(logWeights);
            weights    = weights ./ sum(weights, 2);
            softSym    = weights * constellation;
            confidence = max(weights, [], 2);
        end

        % ESTIMATEOUTPUTNOISE 基于最近星座距离估计输出噪声方差。
        function outputNoiseVar = estimateOutputNoise(rawEst, scaledConstellation, noisePowerLin, outputNoiseKeepRatio)
            numSymbols    = length(rawEst);
            distToNearest = min(abs(rawEst - scaledConstellation.').^2, [], 2);
            sortedDist     = sort(distToNearest, 'ascend');
            numKeep        = max(floor(numSymbols * outputNoiseKeepRatio), 1);
            outputNoiseVar = max(mean(sortedDist(1:numKeep)), noisePowerLin);
        end

        % COMPUTEFEEDBACKSYMBOLS 根据反馈模式生成软/硬反馈符号。
        function [feedbackSym, confidence] = computeFeedbackSymbols( ...
                rawEst, constellation, noiseVar, feedbackMode)
            if strcmpi(feedbackMode, 'hard')
                [~, idx] = min(abs(rawEst - constellation.').^2, [], 2);
                feedbackSym = constellation(idx);
                confidence = ones(size(feedbackSym));
            else
                [feedbackSym, confidence] = GiFreeReceiver.computeSoftSymbols( ...
                    rawEst, constellation, noiseVar);
            end
        end

        % COMPUTECLEANINGSYMBOLS 生成用于导频清洗的数据重构符号。
        function cleaningSymbols = computeCleaningSymbols( ...
                rawEst, constellation, noiseVar, cleaningFeedbackMode)
            if strcmpi(cleaningFeedbackMode, 'hard')
                [~, idx] = min(abs(rawEst - constellation.').^2, [], 2);
                cleaningSymbols = constellation(idx);
            else
                [cleaningSymbols, ~] = GiFreeReceiver.computeSoftSymbols( ...
                    rawEst, constellation, noiseVar);
            end
        end

        % SPLITPHASEBUDGET 按总迭代轮次分配粗锁定/迭代精炼/最终定型阶段预算。
        function [numPhase1, numPhase2, numPostDecision] = splitPhaseBudget( ...
                totalIter, enableDfp, prioritizePhase1Budget, ...
                phase1BudgetRatio, phase1BudgetRatioAlt)
            if nargin < 3
                prioritizePhase1Budget = false;
            end

            if totalIter <= 1
                numPhase1 = 1;
                numPhase2 = 0;
                numPostDecision = 0;
            elseif totalIter == 2
                numPhase1 = 1;
                numPhase2 = 1;
                numPostDecision = 0;
            else
                numPostDecision = double(enableDfp);
                remainIter = totalIter - numPostDecision;
                if prioritizePhase1Budget && remainIter >= 5
                    numPhase1 = max(round(remainIter * phase1BudgetRatio), 4);
                else
                    numPhase1 = max(round(remainIter * phase1BudgetRatioAlt), 1);
                end
                numPhase1 = min(numPhase1, remainIter - 1);
                numPhase2 = remainIter - numPhase1;
            end
        end

        % PHASEKEEPRATIO 计算当前阶段置信门控保留比例。
        function keepRatio = phaseKeepRatio( ...
                iter, totalIter, phaseName, enableConfidenceGating, cleaningGateMode, confidence, ...
                adaptConfScaleDiv, adaptConfScaleMin, adaptConfScaleMax, adaptConfKeepCap)
            if ~enableConfidenceGating
                keepRatio = 1.0;
                return;
            end

            if totalIter <= 1
                progress = 1.0;
            else
                progress = (iter - 1) / (totalIter - 1);
            end

            switch lower(phaseName)
                case 'phase1'
                    ratioMin = 0.0;
                    ratioMax = 0.40;
                case 'phase2'
                    ratioMin = 0.35;
                    ratioMax = 0.85;
                otherwise
                    ratioMin = 0.85;
                    ratioMax = 1.00;
            end

            keepRatio = ratioMin + (ratioMax - ratioMin) * progress;
            if strcmpi(cleaningGateMode, 'adaptive')
                meanConf = mean(confidence);
                adaptScale = min(max((1.0 - meanConf) / adaptConfScaleDiv, adaptConfScaleMin), adaptConfScaleMax);
                keepRatio = min(max(keepRatio * adaptScale, ratioMin), adaptConfKeepCap);
            end
        end

        % APPLYCONFIDENCEGATE 按置信度保留高可靠符号，抑制不确定分量。
        function gatedSym = applyConfidenceGate(symbols, confidence, keepRatio)
            numSym = length(symbols);
            if numSym == 0
                gatedSym = symbols;
                return;
            end

            if keepRatio <= 0
                gatedSym = zeros(size(symbols));
                return;
            end

            if keepRatio >= 0.999
                gatedSym = symbols;
                return;
            end

            nKeep = max(1, ceil(numSym * keepRatio));
            [~, order] = sort(confidence, 'descend');
            keepIdx = order(1:nKeep);
            gatedSym = zeros(size(symbols));
            gatedSym(keepIdx) = symbols(keepIdx);
        end

        % COMPUTEPHASE1DETECTIONSTRICTNESS 计算 Phase 1 渐进检测强度。
        function effectiveDetectionStrictness = computePhase1DetectionStrictness( ...
                iter, totalIter, enableProgressiveCfar, initScale, finalScale)
            if ~enableProgressiveCfar
                effectiveDetectionStrictness = 1.0;
                return;
            end

            if totalIter <= 1
                progress = 1.0;
            else
                progress = (iter - 1) / (totalIter - 1);
            end

            effectiveDetectionStrictness = initScale + (finalScale - initScale) * progress;
            effectiveDetectionStrictness = max(effectiveDetectionStrictness, 0);
        end

        % INITIALIZEPATHHISTORY 用初始 OMP 结果创建路径稳定度历史。
        function pathHistory = initializePathHistory(estPaths)
            pathHistory = GiFreeReceiver.createEmptyPathHistory();
            if isempty(estPaths)
                return;
            end

            numPaths = size(estPaths, 1);
            pathHistory(1, numPaths) = struct( ...
                'delay', 0, 'doppler', 0, 'hitCount', 0, 'gain', 0);
            for pathIdx = 1:numPaths
                pathHistory(pathIdx).delay = estPaths(pathIdx, 1);
                pathHistory(pathIdx).doppler = estPaths(pathIdx, 2);
                pathHistory(pathIdx).hitCount = 1;
                pathHistory(pathIdx).gain = estPaths(pathIdx, 3);
            end
        end

        % UPDATEPATHHISTORY 用本轮 OMP 结果更新路径稳定度历史。
        function pathHistory = updatePathHistory(pathHistory, estPaths, dopplerTolerance)
            if isempty(pathHistory)
                pathHistory = GiFreeReceiver.initializePathHistory(estPaths);
                return;
            end

            if isempty(estPaths)
                for historyIdx = 1:numel(pathHistory)
                    pathHistory(historyIdx).hitCount = pathHistory(historyIdx).hitCount - 1;
                end
                pathHistory = pathHistory([pathHistory.hitCount] > 0);
                return;
            end

            originalHistoryCount = numel(pathHistory);
            matchedHistory = false(originalHistoryCount, 1);
            for pathIdx = 1:size(estPaths, 1)
                candidatePath = estPaths(pathIdx, :);
                matchedIdx = GiFreeReceiver.findPathHistoryMatch( ...
                    pathHistory, matchedHistory, candidatePath, dopplerTolerance);
                if matchedIdx > 0
                    pathHistory(matchedIdx).doppler = candidatePath(2);
                    pathHistory(matchedIdx).gain = candidatePath(3);
                    pathHistory(matchedIdx).hitCount = pathHistory(matchedIdx).hitCount + 1;
                    matchedHistory(matchedIdx) = true;
                else
                    pathHistory(end + 1) = struct( ... %#ok<AGROW>
                        'delay', candidatePath(1), ...
                        'doppler', candidatePath(2), ...
                        'hitCount', 1, ...
                        'gain', candidatePath(3));
                end
            end

            for historyIdx = 1:originalHistoryCount
                if ~matchedHistory(historyIdx)
                    pathHistory(historyIdx).hitCount = pathHistory(historyIdx).hitCount - 1;
                end
            end
            pathHistory = pathHistory([pathHistory.hitCount] > 0);
        end

        % EXTRACTSTABLEPATHS 提取满足稳定度门限的路径矩阵。
        function stablePaths = extractStablePaths( ...
                pathHistory, stabilityThreshold, enablePathStabilityGate)
            if isempty(pathHistory)
                stablePaths = zeros(0, 3);
                return;
            end

            if enablePathStabilityGate
                hitCountVec = [pathHistory.hitCount];
                keepMask = hitCountVec >= max(round(stabilityThreshold), 1);
            else
                keepMask = true(1, numel(pathHistory));
            end

            keptHistory = pathHistory(keepMask);
            if isempty(keptHistory)
                stablePaths = zeros(0, 3);
                return;
            end

            stablePaths = zeros(numel(keptHistory), 3);
            for pathIdx = 1:numel(keptHistory)
                stablePaths(pathIdx, :) = [ ...
                    keptHistory(pathIdx).delay, ...
                    keptHistory(pathIdx).doppler, ...
                    keptHistory(pathIdx).gain];
            end
        end

        function restoreDirichletRadius(cfg, dirichletRadius)
            cfg.DirichletRadius = dirichletRadius;
        end

        function [prunedPaths, improveRatioVec] = pruneBootstrapPathsByResidual( ...
                pathSet, residualNormVec, maxBootstrapPaths, improveThreshold)
            if isempty(pathSet)
                prunedPaths = pathSet;
                improveRatioVec = zeros(0, 1);
                return;
            end

            improveRatioVec = GiFreeReceiver.computeBootstrapImproveRatios(residualNormVec);
            numDetected = size(pathSet, 1);
            keepCount = min(numDetected, max(maxBootstrapPaths, 1));

            for pathIdx = 2:min(numDetected, numel(improveRatioVec))
                if improveRatioVec(pathIdx) < improveThreshold
                    keepCount = pathIdx - 1;
                    break;
                end
            end

            keepCount = max(min(keepCount, numDetected), 1);
            prunedPaths = pathSet(1:keepCount, :);
        end

        function improveRatioVec = computeBootstrapImproveRatios(residualNormVec)
            if numel(residualNormVec) <= 1
                improveRatioVec = zeros(0, 1);
                return;
            end

            numPaths = numel(residualNormVec) - 1;
            improveRatioVec = zeros(numPaths, 1);
            for pathIdx = 1:numPaths
                prevResidual = residualNormVec(pathIdx);
                nextResidual = residualNormVec(pathIdx + 1);
                improveRatioVec(pathIdx) = max(prevResidual - nextResidual, 0) / ...
                    max(prevResidual, 1e-12);
            end
        end

        % BUILDOBSERVATIONWEIGHTS 构建 OMP 观测加权，降低不确定泄漏影响。
        function obsWeights = buildObservationWeights( ...
                effectiveChannel, numSc, dataPos1, cleaningSymbols, gatedCleaningSymbols, useWeightedPilotMetric)
            if ~useWeightedPilotMetric
                obsWeights = ones(numSc, 1);
                return;
            end

            uncertainSymbols = cleaningSymbols - gatedCleaningSymbols;
            if ~any(abs(uncertainSymbols) > 0)
                obsWeights = ones(numSc, 1);
                return;
            end

            uncertainFrame = zeros(numSc, 1);
            uncertainFrame(dataPos1) = uncertainSymbols;
            leakProxy = abs(effectiveChannel * uncertainFrame).^2;
            scale = max(mean(leakProxy), 1e-12);
            obsWeights = 1 ./ (1 + leakProxy / scale);
            obsWeights = max(obsWeights, 0.15);
        end

        % BUILDTXFRAME 由导频与数据符号组合构造完整估计发射帧。
        function txFrame = buildTxFrame(numSc, pilotPos1, dataPos1, ampPP, pilotSeq, dataSymbols)
            txFrame = zeros(numSc, 1);
            txFrame(pilotPos1) = ampPP * pilotSeq;
            txFrame(dataPos1) = dataSymbols;
        end

        % BUILDDATAFRAME 仅在数据子载波填充符号，导频位置为零。
        function dataFrame = buildDataFrame(numSc, dataPos1, dataSymbols)
            dataFrame = zeros(numSc, 1);
            dataFrame(dataPos1) = dataSymbols;
        end

        % SELECTOUTPUTNOISEVAR 根据开关选择固定或自适应输出噪声方差。
        function outputNoiseVar = selectOutputNoiseVar( ...
                rawEst, scaledConstellation, noisePowerLin, useAdaptiveOutputNoise, outputNoiseKeepRatio)
            if useAdaptiveOutputNoise
                outputNoiseVar = GiFreeReceiver.estimateOutputNoise( ...
                    rawEst, scaledConstellation, noisePowerLin, outputNoiseKeepRatio);
            else
                outputNoiseVar = noisePowerLin;
            end
        end

        % HASSUPPORTCHANGE 判断两次路径估计支撑集合是否变化。
        function changed = hasSupportChange(oldPaths, newPaths)
            if isempty(oldPaths) && isempty(newPaths)
                changed = false;
                return;
            end
            if size(oldPaths, 1) ~= size(newPaths, 1)
                changed = true;
                return;
            end

            oldSig = [round(oldPaths(:, 1)), round(oldPaths(:, 2))];
            newSig = [round(newPaths(:, 1)), round(newPaths(:, 2))];
            oldSig = sortrows(oldSig, [1 2]);
            newSig = sortrows(newSig, [1 2]);
            changed = ~isequal(oldSig, newSig);
        end

        % PRUNEWEAKPATHS 按增益幅度裁剪弱路径，限制最大路径数。
        function prunedPaths = pruneWeakPaths(pathSet, maxPaths)
            if size(pathSet, 1) <= maxPaths
                prunedPaths = pathSet;
                return;
            end

            [~, order] = sort(abs(pathSet(:, 3)), 'descend');
            keepIdx    = sort(order(1:maxPaths), 'ascend');
            prunedPaths = pathSet(keepIdx, :);
        end

        % CREATEEMPTYPATHHISTORY 构造空的路径稳定度历史结构体。
        function pathHistory = createEmptyPathHistory()
            pathHistory = struct('delay', {}, 'doppler', {}, 'hitCount', {}, 'gain', {});
        end

        % FINDPATHHISTORYMATCH 按 delay 与 Doppler 容差寻找历史匹配项。
        function matchedIdx = findPathHistoryMatch( ...
                pathHistory, matchedHistory, candidatePath, dopplerTolerance)
            matchedIdx = 0;
            bestDopplerDiff = inf;
            bestHitCount = -inf;
            roundedCandidateDoppler = round(candidatePath(2));

            for historyIdx = 1:numel(matchedHistory)
                if matchedHistory(historyIdx)
                    continue;
                end
                if pathHistory(historyIdx).delay ~= candidatePath(1)
                    continue;
                end

                dopplerDiff = abs(round(pathHistory(historyIdx).doppler) - roundedCandidateDoppler);
                if dopplerDiff > dopplerTolerance
                    continue;
                end

                if dopplerDiff < bestDopplerDiff || ...
                        (dopplerDiff == bestDopplerDiff && ...
                        pathHistory(historyIdx).hitCount > bestHitCount)
                    matchedIdx = historyIdx;
                    bestDopplerDiff = dopplerDiff;
                    bestHitCount = pathHistory(historyIdx).hitCount;
                end
            end
        end

    end

end


