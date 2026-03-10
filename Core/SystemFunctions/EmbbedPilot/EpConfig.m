classdef EpConfig < handle
    % EpConfig: EP-AFDM 系统参数配置 (单导频版)
    %
    % 设计理念:
    %   以 TotalSubcarriers (空口物理帧长) 为顶层输入, 确保与 GI-Free 公平对比.
    %   推导链:
    %     TotalSubcarriers → PrefixLength = MaxPathDelays
    %                      → NumDataSubcarriers = TotalSubcarriers - PrefixLength
    %                      → Chirp, ZP, 导频, 资源映射
    %
    % 导频方案: 仅支持单导频, DAFT 域索引 = ZP + 1

    properties (Access = public)
        BitsPerSymbol      (1, 1) double = 2
        ModulationOrder    (1, 1) double = 4
        WaveformType       (1, 1) string = "AFDM"

        TotalSubcarriers   (1, 1) double = 256
        MaxPathDelays      (1, 1) double = 2
        MaxNormDoppler     (1, 1) double = 2
        DopplerGuard       (1, 1) double = 4
        NumPaths           (1, 1) double = 3

        PilotSnr           (1, 1) double = 35

        ManualZeroPaddingLength (1, 1) double = 0
    end

    properties (SetAccess = private)
        NumDataSubcarriers
        PrefixLength
        ChirpParam1
        ChirpParam2
        ZeroPaddingLength
        TheoreticalZeroPaddingLength
        PilotSequenceLength
        PilotIndex
        PilotSequence
        NumActiveCarriers
        ActiveIndices
    end

    methods

        function obj = EpConfig()
            obj.updateDerivedParams();
        end

        function set.TotalSubcarriers(obj, val)
            obj.TotalSubcarriers = val; obj.updateDerivedParams();
        end

        function set.MaxPathDelays(obj, val)
            obj.MaxPathDelays = val; obj.updateDerivedParams();
        end

        function set.ManualZeroPaddingLength(obj, val)
            obj.ManualZeroPaddingLength = val; obj.updateDerivedParams();
        end

        function updateDerivedParams(obj)
            % 初始化守卫: 属性默认值按声明顺序赋值时, 后续属性仍为 []
            if isempty(obj.TotalSubcarriers) || isempty(obj.MaxPathDelays) || ...
               isempty(obj.MaxNormDoppler)   || isempty(obj.DopplerGuard)
                return;
            end

            obj.PrefixLength = obj.MaxPathDelays;
            N = obj.TotalSubcarriers - obj.PrefixLength;
            obj.NumDataSubcarriers = N;

            if N <= 0
                error('EpConfig:BadN', ...
                    'TotalSubcarriers (%d) 须大于 PrefixLength (%d)', ...
                    obj.TotalSubcarriers, obj.PrefixLength);
            end

            obj.ChirpParam1 = (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) / (2 * N);
            obj.ChirpParam2 = 1 / (N ^ 2 * 2 * pi);

            if (2 * (obj.MaxNormDoppler + obj.DopplerGuard) * (obj.MaxPathDelays + 1)) ...
                    + obj.MaxPathDelays > N
                error('EpConfig:Orth', '子载波不满足正交条件');
            end

            zpTheory = ceil((obj.MaxPathDelays + 1) * ...
                (2 * (obj.MaxNormDoppler + obj.DopplerGuard) + 1) - 1);
            obj.TheoreticalZeroPaddingLength = zpTheory;

            if obj.ManualZeroPaddingLength > 0
                obj.ZeroPaddingLength = obj.ManualZeroPaddingLength;
            else
                obj.ZeroPaddingLength = zpTheory;
            end

            % 单导频
            obj.PilotSequenceLength = 1;
            obj.PilotIndex = obj.ZeroPaddingLength + 1;
            obj.PilotSequence = 1;

            pilotSpan = 2 * obj.ZeroPaddingLength + obj.PilotSequenceLength;
            obj.NumActiveCarriers = N - pilotSpan;
            obj.ActiveIndices = (pilotSpan + 1):N;

            if obj.NumActiveCarriers <= 0
                error('EpConfig:NoData', '开销 (%d) 超过子载波 (%d)', pilotSpan, N);
            end
        end

    end

end