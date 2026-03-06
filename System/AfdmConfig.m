classdef AfdmConfig < handle
    % AfdmConfig: AFDM 系统全局参数配置与资源映射类

    properties (Access = public)
        % --- 基础通信与波形参数 ---
        BitsPerSymbol (1, 1) double = 2
        ModulationOrder (1, 1) double = 4
        WaveformType (1, 1) string = "AFDM"

        % --- 物理层帧结构与信道参数 ---
        NumDataSubcarriers (1, 1) double = 256
        MaxNormDoppler (1, 1) double = 2
        MaxPathDelays (1, 1) double = 5
        DopplerGuard (1, 1) double = 4
        NumPaths (1, 1) double = 3

        % --- 导频配置 ---
        PilotType (1, 1) string = "CAZAC" % "Single" 或 "CAZAC"
        PilotSnr (1, 1) double = 35 % 导频信噪比 (dB)

        PulseShapingType = "Rectangular" % 默认矩形窗
    end

    % 由内部自动计算更新
    properties (SetAccess = private)
        % --- 衍生参数 ---
        PrefixLength
        TotalSubcarriers
        ChirpParam1
        ChirpParam2
        ZeroPaddingLength

        % --- 资源映射参数 ---
        PilotSequenceLength
        PilotIndex
        PilotSequence
        NumActiveCarriers
        ActiveIndices

        % --- 窗函数参数 ---
        PulseShapingWindow % 储存生成的列向量窗函数
        TukeyRollOff = 0.2
    end

    methods
        % --- 构造函数 ---
        function obj = AfdmConfig()
            obj.updateDerivedParams();
        end

        % --- 属性监听与自动更新 ---
        % 当用户修改核心参数时，自动触发衍生参数的重新计算
        function set.NumDataSubcarriers(obj, val)
            obj.NumDataSubcarriers = val;
            obj.updateDerivedParams();
        end

        function set.PilotType(obj, val)
            obj.PilotType = val;
            obj.updateDerivedParams();
        end

        % --- 核心计算逻辑 ---
        function updateDerivedParams(obj)
            % 基础物理参数
            obj.PrefixLength = obj.MaxPathDelays;
            obj.TotalSubcarriers = obj.NumDataSubcarriers + obj.PrefixLength;

            % AFDM 核心参数 (c1, c2)
            obj.ChirpParam1 = (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) / (2 * obj.NumDataSubcarriers);
            obj.ChirpParam2 = 1 / (obj.NumDataSubcarriers ^ 2 * 2 * pi);

            % 正交性校验
            if (2 * (obj.MaxNormDoppler + obj.DopplerGuard) * (obj.MaxPathDelays + 1)) + obj.MaxPathDelays > obj.NumDataSubcarriers
                error('AfdmConfig:OrthogonalityViolated', '子载波不满足正交条件，请增大 NumDataSubcarriers 或减小延时多普勒参数。');
            end

            % ZP 长度计算
            zp = ((obj.MaxPathDelays + 1) * (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) - 1); % 这里缩减ZP只是为了测试CAZAC序列的稳定性
            obj.ZeroPaddingLength = ceil(zp);

            % 导频参数与序列初始化
            obj.initPilotSequence();

            % 有效数据映射
            pilotSpan = 2 * obj.ZeroPaddingLength + obj.PilotSequenceLength;
            obj.NumActiveCarriers = obj.NumDataSubcarriers - pilotSpan;
            obj.ActiveIndices = (pilotSpan + 1):obj.NumDataSubcarriers;

            if obj.NumActiveCarriers <= 0
                error('AfdmConfig:InsufficientCarriers', '导频和保护间隔开销太大(总计%d)，导致有效子载波数不够用!', pilotSpan);
            end

            % 更新窗函数
            obj.updatePulseShapingWindow();
        end

    end

    methods (Access = private)
        % --- 导频序列与位置初始化 ---
        function initPilotSequence(obj)
            obj.PilotIndex = obj.ZeroPaddingLength + 1;

            if upper(obj.PilotType) == "SINGLE"
                obj.PilotSequenceLength = 1;
                obj.PilotSequence = 1 + 1j;

            elseif upper(obj.PilotType) == "CAZAC"
                obj.PilotSequenceLength = 15; % 后续可提取为外部可配参数

                % 预生成基础 CAZAC 序列
                % 若后续进行基于不同根索引的 PAPR 优化对比，可将根索引(此处为1)作为 Config 属性传入
                baseSeq = zadoffChuSeq(1, obj.PilotSequenceLength);

                % 这里仅保存归一化的基础序列，具体的功率缩放(scalePilotPower)留给 Tx/Rx 环节结合噪声功率处理
                obj.PilotSequence = baseSeq;
            else
                error('AfdmConfig:UnknownPilot', '未知的 PilotType: %s', obj.PilotType);
            end

        end

        function updatePulseShapingWindow(obj)
            numSubcarriers = obj.NumDataSubcarriers;

            switch upper(obj.PulseShapingType)
                case "RECTANGULAR"
                    obj.PulseShapingWindow = ones(numSubcarriers, 1);
                case "HAMMING"
                    obj.PulseShapingWindow = hamming(numSubcarriers);
                case "TUKEY"
                    obj.PulseShapingWindow = tukeywin(numSubcarriers, obj.TukeyRollOff);
                    % case "CHEBYSHEV"
                    %     obj.PulseShapingWindow = chebwin(numSubcarriers, obj.ChebyshevSidelobeAttenuation);
                otherwise
                    error('Unknown PulseShapingType');
            end

            % 功率归一化，确保加窗不改变总发射功率
            obj.PulseShapingWindow = obj.PulseShapingWindow / norm(obj.PulseShapingWindow) * sqrt(numSubcarriers);
        end

    end

end
