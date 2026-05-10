function resultStruct = Exp4_FractionalDoppler( ...
        numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
        maxDopplerIdx, dopplerGuard, numPaths, pilotSnrEp, pilotSnrGf, ...
        snrVecDb, figuresDir, trialCountOverride)
% EXP4_FRACTIONALDOPPLER 分数 Doppler 下 EP(xi_nu=1/2/4) 与 ProGuard-CLIP 的对比。
%
%   描述:
%   在线运行 Monte Carlo 仿真, 对同一组随机路径参数并行评测:
%     - EP 三组: xi_nu ∈ {1, 2, 4}, 每组单独配置 EpConfig;
%     - GI-Free ProGuard-CLIP: 单帧 CLIP 迭代接收机。
%   每个 trial 采样一次 (delay, doppler, gain) 后, EP 与 GI-Free 共享该信道参数,
%   分别生成自己的物理/有效信道矩阵, 从而显著降低 BER 估计方差。
%
%   修复要点 (2026-04-24):
%     P0 — pilot Es/N0 对等化: 将 ProGuard-CLIP 的 GiFreeConfig.UseDynamicPilot
%          置为 true, 使 pilot_power = dataSnrLin * 10^(DynamicPilotBaseDb/10)。
%          此时 GI-Free 与 EP 的 pilot/noise 均为 baseDb + dataSNR_dB, 消除了
%          原图 (fixed pilot = 35 dB) 在高 SNR 区间因导频-信道估计被钉死而
%          无法咬住 EP 基线的系统性偏差。
%
%   语法:
%   resultStruct = Exp4_FractionalDoppler( ...
%       numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
%       maxDopplerIdx, dopplerGuard, numPaths, pilotSnrEp, pilotSnrGf, ...
%       snrVecDb, figuresDir)
%
%   输入:
%   numSubcarriers  - (double) DAFT 域子载波数; EP 以此作为 NumDataSubcarriers,
%                     GI-Free 以此作为 NumSubcarriers。
%   modulationOrder - (double) 调制阶数 (e.g. 4 → QPSK)。
%   bitsPerSymbol   - (double) 每符号比特数。
%   maxDelaySamples - (double) 最大时延采样 l_max。
%   maxDopplerIdx   - (double) 最大整数 Doppler 索引 alpha_max。
%   dopplerGuard    - (double) GI-Free 的 DopplerGuard (影响 C1/LocStep)。
%   numPaths        - (double) 真实路径数 P。
%   pilotSnrEp      - (double) EP 的 PilotSnrDb (绝对 pilot power dB)。
%   pilotSnrGf      - (double) GI-Free 的 PilotSnrDb, 同时用作 DynamicPilotBaseDb。
%   snrVecDb        - (Nx1)   数据 SNR 扫描点 (dB)。
%   figuresDir      - (char)  图/MAT 输出目录。
%   trialCountOverride - (可选, double 空/标量/向量) trial 数覆盖:
%                        []        - 使用默认自适应调度 (buildExp4TrialCountVec)
%                        标量 K    - 所有 SNR 点均用 K 次 trial (轻量化小跑)
%                        (Nx1) 向量 - 逐 SNR 点覆盖, 必须与 snrVecDb 同长度
%
%   输出:
%   resultStruct 字段:
%       snrVecDb        (Nx1)   扫描 SNR 点
%       xiNuVec         (1x3)   [1 2 4]
%       berEp           (Nx3)   三组 EP BER
%       berAblation     (Nx1)   ProGuard-CLIP BER
%       nmseAblation    (Nx1)   ProGuard-CLIP NMSE (线性)
%       ablationLabels  (1x1)   {'ProGuard-CLIP'}
%       trialCountVec   (Nx1)   每个 SNR 使用的 trial 数
%       elapsedSec      (1,1)   总耗时
%
%   版本历史:
%   2026-04-24 - Aiden - 修复 pilot Es/N0 对等性，恢复在线 MC 执行。

    arguments
        numSubcarriers  (1,1) double {mustBePositive, mustBeInteger}
        modulationOrder (1,1) double {mustBePositive, mustBeInteger}
        bitsPerSymbol   (1,1) double {mustBePositive, mustBeInteger}
        maxDelaySamples (1,1) double {mustBeNonnegative, mustBeInteger}
        maxDopplerIdx   (1,1) double {mustBeNonnegative, mustBeInteger}
        dopplerGuard    (1,1) double {mustBeNonnegative, mustBeInteger}
        numPaths        (1,1) double {mustBePositive, mustBeInteger}
        pilotSnrEp      (1,1) double {mustBeFinite}
        pilotSnrGf      (1,1) double {mustBeFinite}
        snrVecDb        (:,1) double {mustBeFinite}
        figuresDir      (1,:) char
        trialCountOverride      double {mustBeNonnegative} = []
    end

    if ~isfolder(figuresDir)
        mkdir(figuresDir);
    end

    xiNuVec        = [1, 2, 4];
    ablationLabels = {'ProGuard-CLIP'};

    % ---- 构建系统对象 (一次性, 复用) ----
    epSystems = buildEpSystemTrio( ...
        xiNuVec, numSubcarriers, modulationOrder, bitsPerSymbol, ...
        maxDelaySamples, maxDopplerIdx, numPaths, pilotSnrEp);

    gfSystem = GiFreeSystem(buildProGuardConfig( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, ...
        dopplerGuard, numPaths, pilotSnrGf));

    % ---- 预分配结果 ----
    numSnr         = numel(snrVecDb);
    numAblations   = numel(ablationLabels);
    berEp          = zeros(numSnr, numel(xiNuVec));
    berAblation    = zeros(numSnr, numAblations);
    nmseAblation   = zeros(numSnr, numAblations);
    errorCount     = zeros(numSnr, numAblations);
    totalBitsArr   = zeros(numSnr, numAblations);
    trialCountArr  = zeros(numSnr, numAblations);
    ci95Low        = zeros(numSnr, numAblations);
    ci95High       = zeros(numSnr, numAblations);
    trialCountVec  = resolveExp4TrialCountVec(snrVecDb, trialCountOverride);

    fprintf('Exp4_FractionalDoppler: EP (xi_nu = 1/2/4) vs ProGuard-CLIP (Dynamic Pilot, base=%.0f dB)\n', pilotSnrGf);
    fprintf('%-5s %-10s %-10s %-10s %-10s %-11s %-7s\n', ...
        'SNR', 'EP_xi1', 'EP_xi2', 'EP_xi4', 'ProGuard', 'NMSE_PG(dB)', 'Trials');
    totalTic = tic;

    for snrIdx = 1:numSnr
        snrDb     = snrVecDb(snrIdx);
        snrLin    = 10 ^ (snrDb / 10);
        numTrials = trialCountVec(snrIdx);

        epErr  = zeros(numel(xiNuVec), 1);
        epBits = zeros(numel(xiNuVec), 1);
        pgErr  = 0;
        pgBits = 0;
        pgNmseSum = 0;

        for trialIdx = 1:numTrials
            % 共享物理信道参数 (分数 Doppler)
            [pathDelays, pathDopplers, pathGains] = generateRandomChannel( ...
                numPaths, maxDelaySamples, maxDopplerIdx, true);

            % ---- EP 三组 ----
            for epIdx = 1:numel(xiNuVec)
                rEp = epSystems{epIdx}.runTrial( ...
                    snrDb, pathDelays, pathDopplers, pathGains);
                epErr(epIdx)  = epErr(epIdx)  + rEp.bitErrors;
                epBits(epIdx) = epBits(epIdx) + rEp.totalBits;
            end

            % ---- GI-Free ProGuard-CLIP ----
            % 动态导频已在 GiFreeSystem.runTrialWithChannel 内部自动同步:
            %   cfg.CurrentDataSnrLin = dataSnrLin; (参见 GiFreeSystem.m:100-102)
            rGf = gfSystem.runTrialWithChannel( ...
                snrLin, 1.0, pathDelays, pathDopplers, pathGains);
            pgErr     = pgErr     + rGf.bitErrorsSys;
            pgBits    = pgBits    + rGf.totalBits;
            pgNmseSum = pgNmseSum + rGf.mseSystem;
        end

        for epIdx = 1:numel(xiNuVec)
            berEp(snrIdx, epIdx) = berFloor(epErr(epIdx), epBits(epIdx));
        end
        berAblation(snrIdx, 1)  = berFloor(pgErr, pgBits);
        nmseAblation(snrIdx, 1) = pgNmseSum / max(numTrials, 1);

        % 记录原始误码计数与 95% 置信区间（用于后续诊断与 M1 验收）
        errorCount(snrIdx, 1)    = pgErr;
        totalBitsArr(snrIdx, 1)  = pgBits;
        trialCountArr(snrIdx, 1) = numTrials;
        [ci95Low(snrIdx, 1), ci95High(snrIdx, 1)] = computeConfidenceInterval( ...
            pgErr, max(pgBits, 1), 0.05);

        fprintf(' %+4.0f %.2e   %.2e   %.2e   %.2e   %+7.2f    %5d\n', ...
            snrDb, berEp(snrIdx, 1), berEp(snrIdx, 2), berEp(snrIdx, 3), ...
            berAblation(snrIdx, 1), ...
            10 * log10(nmseAblation(snrIdx, 1) + 1e-20), numTrials);
    end

    elapsedSec = toc(totalTic);
    fprintf('Exp4 total elapsed: %.1f s\n', elapsedSec);

    resultStruct = struct( ...
        'snrVecDb',       snrVecDb(:), ...
        'berEp',          berEp, ...
        'berAblation',    berAblation, ...
        'nmseAblation',   nmseAblation, ...
        'xiNuVec',        xiNuVec, ...
        'ablationLabels', {ablationLabels}, ...
        'trialCountVec',  trialCountVec(:), ...
        'elapsedSec',     elapsedSec, ...
        'errorCount',     errorCount, ...
        'totalBits',      totalBitsArr, ...
        'trialCount',     trialCountArr, ...
        'ci95Low',        ci95Low, ...
        'ci95High',       ci95High);

    plotFractionalBenchmark(resultStruct, figuresDir);
    saveResultMat(resultStruct, figuresDir);
end

% ============================================================
%  Helper: EP 三组系统构造 (xi_nu = 1, 2, 4)
% ============================================================

function epSystems = buildEpSystemTrio(xiNuVec, numSubcarriers, modulationOrder, ...
        bitsPerSymbol, maxDelaySamples, maxDopplerIdx, numPaths, pilotSnrDb)
% BUILDEPSYSTEMTRIO 为三种 xi_nu 分别构造一个 EpSystem。
%
%   NOTE: EpConfig 的派生参数由 set.TotalSubcarriers 触发 updateDerivedParams,
%         因此必须先设置 DopplerGuard/MaxDopplerIdx/MaxDelaySamples,
%         最后再赋 TotalSubcarriers, 才能一次性算出 ChirpParam1/PilotIndex/
%         NumActiveCarriers 等。

    epSystems = cell(1, numel(xiNuVec));
    for idx = 1:numel(xiNuVec)
        cfg                        = EpConfig();
        cfg.ModulationOrder        = modulationOrder;
        cfg.BitsPerSymbol          = bitsPerSymbol;
        cfg.NumPaths               = numPaths;
        cfg.PilotSnrDb             = pilotSnrDb;
        cfg.MaxDopplerIdx          = maxDopplerIdx;
        cfg.DopplerGuard           = xiNuVec(idx);
        cfg.MaxDelaySamples        = maxDelaySamples;
        cfg.ManualZeroPaddingLength = 0;
        % 最后赋值, 触发 updateDerivedParams
        cfg.TotalSubcarriers       = numSubcarriers + maxDelaySamples;
        epSystems{idx} = EpSystem(cfg);
    end
end

% ============================================================
%  Helper: ProGuard-CLIP 配置 (含 P0 修复)
% ============================================================

function cfg = buildProGuardConfig(numSubcarriers, modulationOrder, ...
        maxDelaySamples, maxDopplerIdx, dopplerGuard, numPaths, pilotSnrGfDb)
% BUILDPROGUARDCONFIG 构造 GI-Free ProGuard-CLIP 的配置, 动态导频默认打开。

    cfg                    = GiFreeConfig();
    cfg.NumSubcarriers     = numSubcarriers;
    cfg.ModulationOrder    = modulationOrder;
    cfg.MaxDelaySamples    = maxDelaySamples;
    cfg.MaxDopplerIdx      = maxDopplerIdx;
    cfg.NumPaths           = numPaths;
    cfg.NumPathsUpper      = max(numPaths + 3, 6);
    cfg.DopplerGuard       = dopplerGuard;
    cfg.DirichletRadius    = 4;
    cfg.UseFractionalDoppler = true;
    cfg.PilotSnrDb         = pilotSnrGfDb;
    cfg.MaxSicIterations   = 10;

    % ---- P0 修复: 动态导频, pilot Es/N0 ~ dataSnrDb + baseDb ----
    cfg.UseDynamicPilot    = true;
    cfg.DynamicPilotBaseDb = pilotSnrGfDb;

    % ---- CLIP 组件开关 (与 createClipPilotConfigPair 保持一致) ----
    cfg.EnableProgressiveCfar   = true;
    cfg.EnablePathStabilityGate = true;
    cfg.UseWeightedPilotMetric  = true;
    cfg.UseMatchedPilotMetric   = true;
    cfg.EnableDfp               = true;
    cfg.EnableResidualInjection = true;
    cfg.EnableConfidenceGating  = true;
end

% ============================================================
%  Helper: trial 数调度
% ============================================================

function trialVec = resolveExp4TrialCountVec(snrVecDb, override)
% RESOLVEEXP4TRIALCOUNTVEC 三态分发: []→默认自适应, 标量→广播, 向量→逐点覆盖。
%   返回值恒为 numel(snrVecDb)x1 double, 调用端可直接索引。

    numSnr = numel(snrVecDb);
    if isempty(override)
        trialVec = buildExp4TrialCountVec(snrVecDb);
        return;
    end

    if isscalar(override)
        trialVec = round(override) * ones(numSnr, 1);
    elseif numel(override) == numSnr
        trialVec = round(override(:));
    else
        error('Exp4_FractionalDoppler:BadTrialCountOverride', ...
            'trialCountOverride 必须为空/标量/与 snrVecDb 同长度的向量 (当前 %d vs %d)。', ...
            numel(override), numSnr);
    end

    if any(trialVec <= 0)
        error('Exp4_FractionalDoppler:NonPositiveTrialCount', ...
            'trialCountOverride 解析结果含有非正值。');
    end
end

function trialVec = buildExp4TrialCountVec(snrVecDb)
% BUILDEXP4TRIALCOUNTVEC 根据 SNR 分级调整 MC 次数, 保证 BER <= 1e-4 区间的统计分辨率。
%
%   每 trial 产生 (NumSubcarriers-1) * bitsPerSymbol 比特, 对 N=512/QPSK 约 1022 bit;
%   因此 20 dB 目标 BER ~ 1e-4 需要至少 ~1000 trials 才能见到 O(100) 个错码。

    trialVec = 300 * ones(numel(snrVecDb), 1);
    trialVec(snrVecDb <= 5)  = 400;
    trialVec(snrVecDb == 10) = 600;
    trialVec(snrVecDb == 15) = 1200;
    trialVec(snrVecDb >= 20) = 2000;
end

% ============================================================
%  Helper: 结果落盘
% ============================================================

function saveResultMat(resultStruct, figuresDir)
% SAVERESULTMAT 保存结果至 MAT，字段名为 part4Result。

    part4Result = resultStruct; %#ok<NASGU>
    matPath = fullfile(figuresDir, 'Fig4_FractionalDoppler.mat');
    save(matPath, 'part4Result');
    fprintf('  Exp4 结果已保存: %s\n', matPath);
end

% ============================================================
%  Plotter: BER + NMSE 双子图
% ============================================================

function plotFractionalBenchmark(resultStruct, figuresDir)
    style = PlotStyle();
    figHandle = style.newFigure(1080, 430);
    berAxis  = subplot(1, 2, 1, 'Parent', figHandle);
    nmseAxis = subplot(1, 2, 2, 'Parent', figHandle);

    plotBerPanel(berAxis,  resultStruct, style);
    plotNmsePanel(nmseAxis, resultStruct, style);

    sgtitle(figHandle, 'Fractional Doppler: EP and GI-Free ProGuard-CLIP', ...
        'FontName', style.FontName, 'FontSize', style.TitleSize, ...
        'Interpreter', 'none');

    figPath = fullfile(figuresDir, 'Fig4_FractionalDoppler.fig');
    pngPath = fullfile(figuresDir, 'Fig4_FractionalDoppler.png');
    savefig(figHandle, figPath);
    try
        exportgraphics(figHandle, pngPath, 'Resolution', 300);
    catch exception
        warning('Exp4_FractionalDoppler:ExportFailed', ...
            'exportgraphics 失败: %s', exception.message);
    end
end

function plotBerPanel(axHandle, resultStruct, style)
    hold(axHandle, 'on');
    epColors = {style.BlueLight, style.Blue, style.BlueDark};
    epStyles = {'--o', '-.s', '-d'};

    for epIdx = 1:numel(resultStruct.xiNuVec)
        semilogy(axHandle, resultStruct.snrVecDb, ...
            max(resultStruct.berEp(:, epIdx), 1e-6), epStyles{epIdx}, ...
            'Color', epColors{epIdx}, ...
            'LineWidth', 1.6, 'MarkerSize', 6, ...
            'DisplayName', sprintf('EP $\\xi_\\nu=%d$', resultStruct.xiNuVec(epIdx)));
    end

    proIdx = find(strcmp(resultStruct.ablationLabels, 'ProGuard-CLIP'), 1, 'first');
    if isempty(proIdx)
        proIdx = 1;
    end
    semilogy(axHandle, resultStruct.snrVecDb, ...
        max(resultStruct.berAblation(:, proIdx), 1e-6), '-^', ...
        'Color', style.Red, 'LineWidth', 1.9, 'MarkerSize', 7, ...
        'DisplayName', 'GI-Free ProGuard-CLIP');

    style.apply(axHandle, 10);
    set(axHandle, 'YScale', 'log');
    ylim(axHandle, [1e-5, 1]);
    xlabel(axHandle, 'Data SNR (dB)', 'Interpreter', 'latex');
    ylabel(axHandle, 'BER',          'Interpreter', 'latex');
    title(axHandle,  'BER',          'Interpreter', 'latex');
    legendHandle = legend(axHandle, 'Location', 'southwest');
    style.applyLegend(legendHandle, 8);
end

function plotNmsePanel(axHandle, resultStruct, style)
    hold(axHandle, 'on');
    proIdx = find(strcmp(resultStruct.ablationLabels, 'ProGuard-CLIP'), 1, 'first');
    if isempty(proIdx)
        proIdx = 1;
    end

    nmseDb = 10 * log10(resultStruct.nmseAblation(:, proIdx) + 1e-20);
    plot(axHandle, resultStruct.snrVecDb, nmseDb, '-^', ...
        'Color', style.Red, 'LineWidth', 1.9, 'MarkerSize', 7, ...
        'DisplayName', 'GI-Free ProGuard-CLIP');

    style.apply(axHandle, 10);
    xlabel(axHandle, 'Data SNR (dB)',           'Interpreter', 'latex');
    ylabel(axHandle, 'NMSE (dB)',               'Interpreter', 'latex');
    title(axHandle,  'Channel Estimation NMSE', 'Interpreter', 'latex');
    legendHandle = legend(axHandle, 'Location', 'northeast');
    style.applyLegend(legendHandle, 8);
end
