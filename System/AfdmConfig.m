classdef AfdmConfig < handle
    % AfdmConfig: AFDM 系统全局参数配置与资源映射类
    %
    % 整理版 — 保留核心功能:
    %   1. 基础 AFDM 物理参数自动计算 (c1, c2, CPP, ZP)
    %   2. CAZAC / Single 导频序列生成与资源映射
    %   3. ManualZeroPaddingLength: 手动覆盖 ZP (用于频谱效率实验)
    %   4. TheoreticalZeroPaddingLength: 始终保存理论完整值供对比
    %
    % 已移除: 脉冲成形窗 / 帧布局选项 (后续可作为扩展重新引入)

    properties (Access = public)
        % --- 基础通信参数 ---
        BitsPerSymbol (1, 1) double = 2
        ModulationOrder (1, 1) double = 4
        WaveformType (1, 1) string = "AFDM"

        % --- 物理层帧结构 ---
        NumDataSubcarriers (1, 1) double = 256
        MaxNormDoppler (1, 1) double = 2
        MaxPathDelays (1, 1) double = 5
        DopplerGuard (1, 1) double = 4
        NumPaths (1, 1) double = 3

        % --- 导频配置 ---
        PilotType (1, 1) string = "CAZAC" % "Single" 或 "CAZAC"
        PilotSnr (1, 1) double = 35 % 导频信噪比 (dB)

        % --- ZP 手动覆盖 ---
        % 设为正值时覆盖自动计算的 ZP, 设为 0 时恢复自动 (默认)
        ManualZeroPaddingLength (1, 1) double = 0
    end

    properties (SetAccess = private)
        % --- 衍生参数 (自动计算) ---
        PrefixLength
        TotalSubcarriers
        ChirpParam1
        ChirpParam2
        ZeroPaddingLength
        TheoreticalZeroPaddingLength % 始终保存理论完整值

        % --- 资源映射 ---
        PilotSequenceLength
        PilotIndex
        PilotSequence
        NumActiveCarriers
        ActiveIndices
    end

    methods

        function obj = AfdmConfig()
            obj.updateDerivedParams();
        end

        % --- set 监听: 修改关键参数后自动重算 ---
        function set.NumDataSubcarriers(obj, val)
            obj.NumDataSubcarriers = val; obj.updateDerivedParams();
        end

        function set.PilotType(obj, val)
            obj.PilotType = val; obj.updateDerivedParams();
        end

        function set.ManualZeroPaddingLength(obj, val)
            obj.ManualZeroPaddingLength = val; obj.updateDerivedParams();
        end

        function updateDerivedParams(obj)
            N = obj.NumDataSubcarriers;

            % 基础物理参数
            obj.PrefixLength = obj.MaxPathDelays;
            obj.TotalSubcarriers = N + obj.PrefixLength;

            % AFDM chirp 参数
            obj.ChirpParam1 = (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) / (2 * N);
            obj.ChirpParam2 = 1 / (N ^ 2 * 2 * pi);

            % 正交性校验
            if (2 * (obj.MaxNormDoppler + obj.DopplerGuard) * (obj.MaxPathDelays + 1)) + obj.MaxPathDelays > N
                error('AfdmConfig:OrthogonalityViolated', ...
                '子载波不满足正交条件，请增大 NumDataSubcarriers 或减小延时多普勒参数。');
            end

            % ZP 长度: 理论值始终计算, 实际值可手动覆盖
            zpTheory = ceil((obj.MaxPathDelays + 1) * (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) - 1);
            obj.TheoreticalZeroPaddingLength = zpTheory;

            if obj.ManualZeroPaddingLength > 0
                obj.ZeroPaddingLength = obj.ManualZeroPaddingLength;
            else
                obj.ZeroPaddingLength = zpTheory;
            end

            % 导频初始化
            obj.initPilotSequence();

            % 资源映射: [保护带 | 导频 | 保护带 | ═══ 数据 ═══]
            pilotSpan = 2 * obj.ZeroPaddingLength + obj.PilotSequenceLength;
            obj.NumActiveCarriers = N - pilotSpan;
            obj.ActiveIndices = (pilotSpan + 1):N;

            if obj.NumActiveCarriers <= 0
                error('AfdmConfig:InsufficientCarriers', ...
                    '保护带+导频开销 (%d) 导致有效子载波不足!', pilotSpan);
            end

        end

    end

    methods (Access = private)

        function initPilotSequence(obj)
            obj.PilotIndex = obj.ZeroPaddingLength + 1;

            if upper(obj.PilotType) == "SINGLE"
                obj.PilotSequenceLength = 1;
                obj.PilotSequence = 1 + 1j;
            elseif upper(obj.PilotType) == "CAZAC"
                obj.PilotSequenceLength = 15;
                obj.PilotSequence = zadoffChuSeq(1, obj.PilotSequenceLength);
            else
                error('AfdmConfig:UnknownPilot', '未知的 PilotType: %s', obj.PilotType);
            end

        end

    end

end
