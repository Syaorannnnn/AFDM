classdef GiFreeConfig < handle
% GIFREECONFIG - GI-Free 系统配置对象与派生参数集合
%
%   描述:
%   维护 GI-Free 架构的基础参数与依赖参数（Dependent properties），包含
%   DAFT chirp 参数、导频/数据位置、导频功率映射以及兼容旧命名别名等逻辑。
%
%   语法:
%   cfg = GiFreeConfig();
%   cfg.NumSubcarriers = 512;
%   cfg.validate();
%
%   输出:
%   cfg - (GiFreeConfig) 可直接传入 GiFreeSystem / GiFreeTransmitter / GiFreeReceiver。
%
%   NOTE:
%   PilotPos0/DataPos0 为 0-based 语义；PilotPos1/DataPos1 为 MATLAB 1-based 索引。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    properties
        NumSubcarriers
        ModulationOrder
        MaxDelaySamples
        MaxDopplerIdx
        NumPaths
        DopplerGuard
        DirichletRadius      % Dirichlet 核在 Doppler 维的半展宽
        PilotSnrDb
        MaxSicIterations
        NumPathsUpper
        UseFractionalDoppler
        EarlyPathSearchMode = 'greedy' % 'greedy' | 'beam'
        BeamWidth = 2
        BeamDepth = 2
        BeamStrategyMode = 'fastAdaptive' % 'fastAdaptive' | 'performance' | 'fixed'
        PerformanceBeamWidth = 2
        PerformanceBeamDepth = 2
        BeamExpandFactor = 2
        BeamMinImproveRatio = 0.005
        BeamAdaptiveImproveRatio = 0.02
        BeamAdaptiveSpreadThreshold = 4
        BeamUncertainMetricRatio = 1.25
        EnableBeamSignaturePrune = true
        FeedbackMode = 'soft'      % 'soft' | 'hard'
        CleaningFeedbackMode = 'soft' % 'soft' | 'hard'
        CleaningGateMode = 'topk'  % 'topk' | 'adaptive'
        UseAdaptiveOutputNoise = true
        EnableDfp = true
        EnableResidualInjection = true
        EnableConfidenceGating = true
        % Phase 0 bootstrap 已退役，以下参数仅为兼容旧脚本保留。
        EnablePhase0Bootstrap = false
        EnablePhase0LocalRefine = false
        Phase0CleanRatio = 0.08
        BootstrapImproveRatioThreshold = 0.03
        MaxBootstrapPaths = 5
        EnableProgressiveCfar = true
        ProgressiveCfarInitScale = 3.0
        ProgressiveCfarFinalScale = 1.0
        EnablePathStabilityGate = true
        PathStabilityThreshold = 2
        PathStabilityDopplerTolerance = 1
        UseWeightedPilotMetric = true
        UseMatchedPilotMetric = false % 使用完整复合响应归一化匹配作为 OMP 粗搜统计量
        EnableCfarDiagnostic = false  % 输出 metric/threshold 诊断样本，不改变检测流程
        PilotClusterWidth = 2
        % 动态导频模式：导频功率跟踪数据 SNR，使 pilot/noise 恒定（与 EP 行为对齐）
        UseDynamicPilot = false
        DynamicPilotBaseDb = 35   % 动态模式下固定的 pilot/noise 比 (dB)
        CurrentDataSnrLin = 1     % 当前 trial 的数据 SNR 线性值（由 GiFreeSystem 外部设置）
        CfarFalseAlarmProb = 1e-3 % CFAR 目标虚警率
        EnableCfarPilotPowerNormalization = false
        CfarPilotPowerCapDb = 45   % CFAR 门限参考导频功率上限 (dB)
        % OMP 正则化参数
        OmpRegUpperBound = 0.15       % Phase0/1/2 OMP 正则化上限
        Phase1OmpRegUpperBound = 0.45 % Phase1 OMP 正则化上限（渐进 CFAR 时）
        ResidualPathThresholdScale = 0.85 % 残差路径检测门限缩放
        UnionImproveRatio = 0.995     % 合并路径改善接受比例
        % Phase2 迭代控制
        Phase2AcceptTol = 0.01        % Phase2 接受容限（无支撑集变化）
        Phase2AcceptTolWithChange = 0.02 % Phase2 接受容限（有支撑集变化）
        Phase2ImproveTol = 0.005      % Phase2 改善判定门限
        % 噪声估计截尾
        TrimAlphaInit = 0.15          % 截尾 alpha 初值
        TrimAlphaDecay = 0.02         % 截尾 alpha 每迭代衰减量
        TrimAlphaMin = 0.05           % 截尾 alpha 下限
        % 输出噪声与相位预算
        OutputNoiseKeepRatio = 0.9    % 输出噪声估计保留比例
        Phase1BudgetRatio = 0.40      % Phase1 预算比例（prioritize 模式）
        Phase1BudgetRatioAlt = 0.30   % Phase1 预算比例（普通模式）
        % 自适应置信度缩放
        AdaptConfScaleDiv = 0.35      % 自适应置信度缩放除数
        AdaptConfScaleMin = 0.75      % 自适应置信度缩放下限
        AdaptConfScaleMax = 1.35      % 自适应置信度缩放上限
        AdaptConfKeepCap = 0.95       % 自适应置信度保留比例上限
    end

    properties (Dependent)
        ChirpParam1
        ChirpParam2
        PilotAmplitude
        LocStep
        PilotZoneSpan
        PilotPos0
        DataPos0
        PilotPos1
        DataPos1
        NumDataSymbols
        PilotSequence
        PerPilotAmplitude
        SpreadWidth

        MaxDelaySpread
        MaxDopplerIndex
        PilotPositions
        DataPositions
    end

    methods

        function set.MaxDelaySpread(obj, val)
            GiFreeConfig.warnDeprecatedAlias('MaxDelaySpread', 'MaxDelaySamples');
            obj.MaxDelaySamples = val;
        end

        function val = get.MaxDelaySpread(obj)
            val = obj.MaxDelaySamples;
        end

        function set.MaxDopplerIndex(obj, val)
            GiFreeConfig.warnDeprecatedAlias('MaxDopplerIndex', 'MaxDopplerIdx');
            obj.MaxDopplerIdx = val;
        end

        function val = get.MaxDopplerIndex(obj)
            val = obj.MaxDopplerIdx;
        end

        function val = get.PilotPositions(obj)
            GiFreeConfig.warnDeprecatedAlias('PilotPositions', 'PilotPos0');
            val = obj.PilotPos0;
        end

        function val = get.DataPositions(obj)
            GiFreeConfig.warnDeprecatedAlias('DataPositions', 'DataPos0');
            val = obj.DataPos0;
        end

        function set.SpreadWidth(obj, val)
            GiFreeConfig.warnDeprecatedAlias('SpreadWidth', 'DirichletRadius');
            obj.DirichletRadius = val;
        end

        function val = get.SpreadWidth(obj)
            GiFreeConfig.warnDeprecatedAlias('SpreadWidth', 'DirichletRadius');
            val = obj.DirichletRadius;
        end

        % VALIDATE 校验 Theorem-1 相关正交约束是否满足。
        function validate(obj)
            lhs = 2 * obj.MaxDopplerIdx * obj.MaxDelaySamples ...
                + 2 * obj.MaxDopplerIdx + obj.MaxDelaySamples;
            if lhs >= obj.NumSubcarriers
                error('GiFreeConfig:Theorem1Violation', ...
                    ['正交条件不满足: 2*alpha_max*l_max + 2*alpha_max + l_max = %d >= N = %d\n' ...
                     '请减小 MaxDopplerIdx (%d) 或 MaxDelaySamples (%d)，或增大 NumSubcarriers (%d).'], ...
                    lhs, obj.NumSubcarriers, ...
                    obj.MaxDopplerIdx, obj.MaxDelaySamples, obj.NumSubcarriers);
            end
        end

        function val = get.ChirpParam1(obj)
            val = (2 * (obj.MaxDopplerIdx + obj.DopplerGuard) + 1) ...
                / (2 * obj.NumSubcarriers);
        end

        function val = get.ChirpParam2(obj)
            val = 1 / (2 * obj.NumSubcarriers^2 * pi);
        end

        function val = get.PilotAmplitude(obj)
            if obj.UseDynamicPilot
                % 动态模式: pilot_power = dataSnrLin * 10^(baseDb/10)
                % 使 pilot/noise = baseDb + dataSNR_dB，与 EP 归一化一致
                val = sqrt(obj.CurrentDataSnrLin * 10^(obj.DynamicPilotBaseDb / 10));
            else
                val = sqrt(10^(obj.PilotSnrDb / 10));
            end
        end

        function val = get.LocStep(obj)
            val = 2 * (obj.MaxDopplerIdx + obj.DopplerGuard) + 1;
        end

        function val = get.PilotZoneSpan(obj)
            val = obj.LocStep * (obj.MaxDelaySamples + 1);
        end

        function val = get.PilotPos0(~)
            val = 0;
        end

        function val = get.DataPos0(obj)
            allPos0 = (0:obj.NumSubcarriers - 1).';
            isPilot = ismember(allPos0, obj.PilotPos0);
            val = allPos0(~isPilot);
        end

        function val = get.PilotPos1(obj)
            val = obj.PilotPos0 + 1;
        end

        function val = get.DataPos1(obj)
            val = obj.DataPos0 + 1;
        end

        function val = get.NumDataSymbols(obj)
            val = obj.NumSubcarriers - 1;
        end

        function val = get.PilotSequence(~)
            val = 1;
        end

        function val = get.PerPilotAmplitude(obj)
            val = obj.PilotAmplitude;
        end

        function newConfig = copy(obj)
        % COPY 创建 GiFreeConfig 对象的深拷贝。
        %
        % 输出:
        %   newConfig - (GiFreeConfig) 与当前对象属性值一致的新配置对象。

            newConfig = GiFreeConfig();
            propertyMetaList = metaclass(obj).PropertyList;
            for propertyIdx = 1:numel(propertyMetaList)
                propertyMeta = propertyMetaList(propertyIdx);
                if iscell(propertyMetaList)
                    propertyMeta = propertyMetaList{propertyIdx};
                end

                if propertyMeta.Dependent || propertyMeta.Constant
                    continue;
                end

                propertyName = propertyMeta.Name;
                if isprop(newConfig, propertyName)
                    newConfig.(propertyName) = obj.(propertyName);
                end
            end
        end

    end

    methods (Static, Access = private)

        % WARNDEPRECATEDALIAS 对旧字段名访问给出一次性弃用告警。
        function warnDeprecatedAlias(oldName, newName)
            persistent warnedMap;
            if isempty(warnedMap)
                warnedMap = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            end

            warnKey = [oldName '->' newName];
            if ~isKey(warnedMap, warnKey)
                warning('GiFreeConfig:DeprecatedAlias', ...
                    'Property "%s" is deprecated. Use "%s" instead.', oldName, newName);
                warnedMap(warnKey) = true;
            end
        end

    end

end

