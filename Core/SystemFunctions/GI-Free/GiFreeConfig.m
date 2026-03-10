classdef GiFreeConfig < handle
    % GiFreeConfig  GI-Free AFDM 系统参数配置 (handle 类)
    %
    %   作为 handle 类, 所有持有该对象引用的模块共享同一份配置;
    %   修改任一属性后, 所有模块立即可见最新值, 无需重新传递.
    %
    %   集中管理子载波数、调制阶数、信道参数、导频功率等所有可调参数,
    %   通过 Dependent 属性自动计算 chirp 系数、导频幅度、网格步长等派生量。

    properties
        NumSubcarriers        % 子载波数 N
        ModulationOrder       % QAM 调制阶数 M
        MaxDelaySpread        % 最大时延扩展 l_max
        MaxDopplerIndex       % 最大多普勒指标 k_max
        NumPaths              % 信道路径数 P
        DopplerGuard          % 多普勒保护间隔 xi_nu
        SpreadWidth           % Dirichlet 扩展半宽 k_nu
        PilotSnrDb            % 导频信噪比 (dB)
        MaxSicIterations      % 总 SIC 迭代次数 (Phase1 + Phase2 + PostDecision)
        NumPathsUpper         % OMP 搜索路径上限
        UseFractionalDoppler  % 是否启用分数多普勒
    end

    properties (Dependent)
        ChirpParam1     % DAFT 线性 chirp 参数 c1
        ChirpParam2     % DAFT 二次 chirp 参数 c2
        PilotAmplitude  % 导频符号幅度 Ap
        LocStep         % 时延-多普勒网格步长
    end

    methods

        function obj = GiFreeConfig()
            % 默认参数 — 针对最终系统优化
            obj.NumSubcarriers       = 256;
            obj.ModulationOrder      = 4;
            obj.MaxDelaySpread       = 2;
            obj.MaxDopplerIndex      = 2;
            obj.NumPaths             = 3;
            obj.DopplerGuard         = 4;
            obj.SpreadWidth          = 4;
            obj.PilotSnrDb           = 45;
            obj.MaxSicIterations     = 10;
            obj.NumPathsUpper        = 4;
            obj.UseFractionalDoppler = true;
        end

        function val = get.ChirpParam1(obj)
            val = (2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1) / (2 * obj.NumSubcarriers);
        end

        function val = get.ChirpParam2(obj)
            val = 1 / (2 * obj.NumSubcarriers^2 * pi);
        end

        function val = get.PilotAmplitude(obj)
            val = sqrt(10^(obj.PilotSnrDb / 10));
        end

        function val = get.LocStep(obj)
            val = 2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1;
        end

    end

end
