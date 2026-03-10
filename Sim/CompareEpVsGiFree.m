%% CompareEpVsGiFree.m
%  EP vs GI-Free AFDM 公平对比仿真 (单导频版)
%
%  公平对比设计:
%    physicalFrameLength = 两系统空口发射的总采样数
%      EP:      发射 physicalFrameLength 个采样 (含 CPP)
%               → N_data = physicalFrameLength - L_cp, 再扣 ZP + 导频
%      GI-Free: 发射 physicalFrameLength 个采样 (无 CP)
%               → 仅扣 1 个导频
%
%  BER = 0 的处理:
%    高 SNR 下 Monte Carlo 可能观测到零错误.
%    绘图时用 berFloor = 0.5 / totalBitsPerSnrPoint 作为置信下界,
%    表示"在此试验量下无法分辨低于此值的 BER".

clear; clc; close all;

%% ================= 全局参数 =================
physicalFrameLength = 512;
modulationOrder     = 4;
bitsPerSymbol       = log2(modulationOrder);
maxDelaySpread      = 2;
maxDopplerIndex     = 2;
numPaths            = 3;
dopplerGuard        = 4;
snrRange            = 0:5:20;
numTrials           = 1e3;

rng(42);

%% ================= 配置 EP =================
fprintf('====== EP-AFDM (Single Pilot) ======\n');
epCfg = EpConfig();
epCfg.TotalSubcarriers = physicalFrameLength;
epCfg.ModulationOrder = modulationOrder;
epCfg.BitsPerSymbol = bitsPerSymbol;
epCfg.MaxPathDelays = maxDelaySpread;
epCfg.MaxNormDoppler = maxDopplerIndex;
epCfg.NumPaths = numPaths;
epCfg.DopplerGuard = dopplerGuard;
epCfg.PilotSnr = 35;

epSys = EpSystem(epCfg);

fprintf('  空口帧长   = %d\n', epCfg.TotalSubcarriers);
fprintf('  CPP        = %d\n', epCfg.PrefixLength);
fprintf('  N_data     = %d\n', epCfg.NumDataSubcarriers);
fprintf('  ZP (×2)    = %d\n', 2 * epCfg.ZeroPaddingLength);
fprintf('  导频       = %d (单导频)\n', epCfg.PilotSequenceLength);
fprintf('  有效数据   = %d\n', epCfg.NumActiveCarriers);

%% ================= 配置 GI-Free =================
fprintf('\n====== GI-Free AFDM (Single Pilot) ======\n');
gfCfg = GiFreeConfig();
gfCfg.NumSubcarriers = physicalFrameLength;
gfCfg.ModulationOrder = modulationOrder;
gfCfg.MaxDelaySpread = maxDelaySpread;
gfCfg.MaxDopplerIndex = maxDopplerIndex;
gfCfg.NumPaths = numPaths;
gfCfg.DopplerGuard = dopplerGuard;
gfCfg.SpreadWidth = 4;
gfCfg.PilotSnrDb = 45;
gfCfg.MaxSicIterations = 10;
gfCfg.NumPathsUpper = 4;
gfCfg.UseFractionalDoppler = true;

gfSys = GiFreeSystem(gfCfg);

fprintf('  空口帧长   = %d\n', gfCfg.NumSubcarriers);
fprintf('  保护间隔   = 0\n');
fprintf('  有效数据   = %d\n', gfCfg.NumSubcarriers - 1);

%% ================= 频谱效率 =================
fprintf('\n====== 频谱效率 ======\n');

epDataSymbols = epCfg.NumActiveCarriers;
gfDataSymbols = gfCfg.NumSubcarriers - 1;

epSe = (epDataSymbols * bitsPerSymbol) / physicalFrameLength;
gfSe = (gfDataSymbols * bitsPerSymbol) / physicalFrameLength;

epOverhead = 1 - epDataSymbols / physicalFrameLength;
gfOverhead = 1 - gfDataSymbols / physicalFrameLength;
seGain = (gfSe - epSe) / epSe * 100;

fprintf('  EP:      %3d / %d → %.4f bps/Hz (开销 %.1f%%)\n', ...
    epDataSymbols, physicalFrameLength, epSe, epOverhead * 100);
fprintf('  GI-Free: %3d / %d → %.4f bps/Hz (开销 %.1f%%)\n', ...
    gfDataSymbols, physicalFrameLength, gfSe, gfOverhead * 100);
fprintf('  提升: +%.1f%%\n', seGain);

%% ================= Monte Carlo =================
fprintf('\n====== Monte Carlo (%d trials/SNR) ======\n', numTrials);

numSnr = length(snrRange);
epBerEst = zeros(numSnr, 1);   epBerPerf = zeros(numSnr, 1);
gfBerEst = zeros(numSnr, 1);   gfBerPerf = zeros(numSnr, 1);
gfNmse = zeros(numSnr, 1);
epTotalBitsPerSnr = 0;
gfTotalBitsPerSnr = 0;

for si = 1:numSnr
    snrDb = snrRange(si);
    snrLin = 10 ^ (snrDb / 10);
    noisePow = 1 / snrLin;

    epErrE = 0; epErrP = 0; epBitsAcc = 0;
    gfErrS = 0; gfErrR = 0; gfBitsAcc = 0;
    gfMseAcc = 0;

    for t = 1:numTrials
        [dl, dp, gn] = generateChannel(maxDelaySpread, maxDopplerIndex, numPaths);

        % --- EP ---
        rEst = epSys.runTrial(snrDb, dl, dp, gn);
        epErrE = epErrE + rEst.bitErrors;
        epBitsAcc = epBitsAcc + rEst.totalBits;

        rPerf = epSys.runTrialPerfectCsi(snrDb, dl, dp, gn);
        epErrP = epErrP + rPerf.bitErrors;

        % --- GI-Free ---
        % GI-Free 的 SNR = dataSnrLinear / noisePower
        % 传入 (snrLin, 1) 使得实际 SNR = snrLin / 1 = snrLin
        % 注意: 不能传 (snrLin, 1/snrLin), 否则实际 SNR = snrLin²
        rGf = gfSys.runTrial(snrLin, 1);
        gfErrS = gfErrS + rGf.bitErrorsSys;
        gfErrR = gfErrR + rGf.bitErrorsRef;
        gfBitsAcc = gfBitsAcc + rGf.totalBits;
        gfMseAcc = gfMseAcc + rGf.mseSystem;
    end

    epBerEst(si)  = epErrE / epBitsAcc;
    epBerPerf(si) = epErrP / epBitsAcc;
    gfBerEst(si)  = gfErrS / gfBitsAcc;
    gfBerPerf(si) = gfErrR / gfBitsAcc;
    gfNmse(si)    = gfMseAcc / numTrials;

    if si == 1
        epTotalBitsPerSnr = epBitsAcc;
        gfTotalBitsPerSnr = gfBitsAcc;
    end

    fprintf('SNR=%2ddB | EP: E=%.2e P=%.2e | GF: E=%.2e P=%.2e\n', ...
        snrDb, epBerEst(si), epBerPerf(si), gfBerEst(si), gfBerPerf(si));
end

%% ================= BER = 0 处理 =================
% 置信下界: 在 totalBits 中观测到 0 错误, 真实 BER 的 95% 置信上界约为 3/totalBits
% 这里用 0.5/totalBits 作为绘图地板, 零错误点用空心标记区分
epBerFloor  = 0.5 / epTotalBitsPerSnr;
gfBerFloor  = 0.5 / gfTotalBitsPerSnr;

epBerEstPlot  = max(epBerEst,  epBerFloor);
epBerPerfPlot = max(epBerPerf, epBerFloor);
gfBerEstPlot  = max(gfBerEst,  gfBerFloor);
gfBerPerfPlot = max(gfBerPerf, gfBerFloor);

%% ================= 吞吐量 =================
epThEst  = epSe * (1 - epBerEst);
epThPerf = epSe * (1 - epBerPerf);
gfThEst  = gfSe * (1 - gfBerEst);
gfThPerf = gfSe * (1 - gfBerPerf);

%% ================= 配色 =================
cEpE  = [0.173 0.373 0.604];   % 深蓝
cEpP  = [0.459 0.639 0.831];   % 浅蓝
cGfE  = [0.769 0.275 0.216];   % 深红
cGfP  = [0.918 0.553 0.490];   % 浅红
bgC   = [0.98 0.98 0.98];
gridC = [0.82 0.82 0.82];
ovhC  = [0.75 0.75 0.75; 0.60 0.60 0.60; 0.85 0.70 0.45; 0.35 0.65 0.45];

%% ================= Figure 1: BER =================
fig1 = figure('Position', [80 80 740 540], 'Color', 'w');
ax1 = axes('Parent', fig1); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

semilogy(ax1, snrRange, epBerEstPlot, '-o', 'Color', cEpE, 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', cEpE, 'DisplayName', 'EP (Estimated CSI)');
semilogy(ax1, snrRange, epBerPerfPlot, '--s', 'Color', cEpP, 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', cEpP, 'DisplayName', 'EP (Perfect CSI)');
semilogy(ax1, snrRange, gfBerEstPlot, '-^', 'Color', cGfE, 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', cGfE, 'DisplayName', 'GI-Free (Estimated CSI)');
semilogy(ax1, snrRange, gfBerPerfPlot, '--d', 'Color', cGfP, 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', cGfP, 'DisplayName', 'GI-Free (Perfect CSI)');

% 标注零错误点 (空心标记)
zeroEpE = snrRange(epBerEst == 0);
if ~isempty(zeroEpE)
    semilogy(ax1, zeroEpE, epBerFloor * ones(size(zeroEpE)), 'o', ...
        'Color', cEpE, 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

zeroGfE = snrRange(gfBerEst == 0);
if ~isempty(zeroGfE)
    semilogy(ax1, zeroGfE, gfBerFloor * ones(size(zeroGfE)), '^', ...
        'Color', cGfE, 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% BER 地板虚线
yline(ax1, epBerFloor, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, ...
    'Label', 'MC floor', 'FontSize', 8, 'LabelHorizontalAlignment', 'left');

xlabel(ax1, 'SNR (dB)', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel(ax1, 'Bit Error Rate', 'FontSize', 13, 'FontName', 'Times New Roman');
title(ax1, sprintf('BER — Single Pilot (Frame = %d)', physicalFrameLength), ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax1, 'Location', 'southwest', 'FontSize', 10, 'FontName', 'Times New Roman');
set(ax1, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.6);
ylim(ax1, [epBerFloor / 3, 1]); xlim(ax1, [snrRange(1) snrRange(end)]);

exportgraphics(fig1, 'fig1_ber.png', 'Resolution', 300);
fprintf('\n已保存: fig1_ber.png\n');

%% ================= Figure 2: 频谱效率 + 资源分解 =================
fig2 = figure('Position', [80 80 740 480], 'Color', 'w');

% 子图 1: 频谱效率
ax2a = subplot(1, 2, 1); hold(ax2a, 'on'); grid(ax2a, 'on'); box(ax2a, 'on');
b1 = bar(ax2a, [epSe; gfSe], 0.55);
b1.FaceColor = 'flat'; b1.CData = [cEpE; cGfE]; b1.EdgeColor = 'none';
set(ax2a, 'XTickLabel', {'EP', 'GI-Free'}, 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel(ax2a, 'bps/Hz', 'FontSize', 12, 'FontName', 'Times New Roman');
title(ax2a, 'Spectral Efficiency', 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(ax2a, 'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.5);
text(ax2a, 1, epSe + 0.03, sprintf('%.3f', epSe), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
text(ax2a, 2, gfSe + 0.03, sprintf('%.3f (+%.0f%%)', gfSe, seGain), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman', 'Color', cGfE);

% 子图 2: 资源构成
ax2b = subplot(1, 2, 2); hold(ax2b, 'on'); grid(ax2b, 'on'); box(ax2b, 'on');
epCpp = epCfg.PrefixLength;
epZp = 2 * epCfg.ZeroPaddingLength;
epPil = epCfg.PilotSequenceLength;
epDat = epCfg.NumActiveCarriers;
gfPil = 1;
gfDat = gfCfg.NumSubcarriers - 1;

bm = [epCpp 0; epZp 0; epPil gfPil; epDat gfDat];
b2 = bar(ax2b, bm', 'stacked');
for k = 1:4, b2(k).FaceColor = ovhC(k, :); b2(k).EdgeColor = 'none'; end
set(ax2b, 'XTickLabel', {'EP', 'GI-Free'}, 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel(ax2b, 'Samples', 'FontSize', 12, 'FontName', 'Times New Roman');
title(ax2b, sprintf('Resource (Total = %d)', physicalFrameLength), ...
    'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax2b, {'CPP','ZP','Pilot','Data'}, 'Location', 'north', 'FontSize', 9, 'FontName', 'Times New Roman');
set(ax2b, 'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.5);
ylim(ax2b, [0 physicalFrameLength * 1.12]);

exportgraphics(fig2, 'fig2_efficiency.png', 'Resolution', 300);
fprintf('已保存: fig2_efficiency.png\n');

%% ================= Figure 3: 吞吐量 =================
fig3 = figure('Position', [80 80 740 540], 'Color', 'w');
ax3 = axes('Parent', fig3); hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');

plot(ax3, snrRange, epThEst, '-o', 'Color', cEpE, 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', cEpE, 'DisplayName', 'EP (Est)');
plot(ax3, snrRange, epThPerf, '--s', 'Color', cEpP, 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', cEpP, 'DisplayName', 'EP (Perf)');
plot(ax3, snrRange, gfThEst, '-^', 'Color', cGfE, 'LineWidth', 2, ...
    'MarkerSize', 7, 'MarkerFaceColor', cGfE, 'DisplayName', 'GI-Free (Est)');
plot(ax3, snrRange, gfThPerf, '--d', 'Color', cGfP, 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', cGfP, 'DisplayName', 'GI-Free (Perf)');

yline(ax3, epSe, ':', 'Color', cEpE, 'LineWidth', 1, ...
    'Label', sprintf('EP: %.3f', epSe), 'FontSize', 9, 'LabelHorizontalAlignment', 'left');
yline(ax3, gfSe, ':', 'Color', cGfE, 'LineWidth', 1, ...
    'Label', sprintf('GF: %.3f', gfSe), 'FontSize', 9, 'LabelHorizontalAlignment', 'left');

xlabel(ax3, 'SNR (dB)', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel(ax3, 'Throughput (bps/Hz)', 'FontSize', 13, 'FontName', 'Times New Roman');
title(ax3, sprintf('Throughput (Frame = %d)', physicalFrameLength), ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax3, 'Location', 'southeast', 'FontSize', 10, 'FontName', 'Times New Roman');
set(ax3, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.6);
xlim(ax3, [snrRange(1) snrRange(end)]);

exportgraphics(fig3, 'fig3_throughput.png', 'Resolution', 300);
fprintf('已保存: fig3_throughput.png\n');

%% ================= Figure 4: 综合仪表板 =================
fig4 = figure('Position', [40 40 980 700], 'Color', 'w');

% (a) BER
ax4a = subplot(2, 2, 1); hold(ax4a, 'on'); grid(ax4a, 'on'); box(ax4a, 'on');
semilogy(ax4a, snrRange, epBerEstPlot, '-o', 'Color', cEpE, 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'MarkerFaceColor', cEpE);
semilogy(ax4a, snrRange, gfBerEstPlot, '-^', 'Color', cGfE, 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'MarkerFaceColor', cGfE);
set(ax4a, 'YScale', 'log');
xlabel(ax4a, 'SNR (dB)', 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel(ax4a, 'BER', 'FontSize', 11, 'FontName', 'Times New Roman');
title(ax4a, '(a) BER (Estimated CSI)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax4a, {'EP','GI-Free'}, 'Location', 'southwest', 'FontSize', 9);
set(ax4a, 'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.5);
ylim(ax4a, [epBerFloor / 3, 1]);

% (b) 吞吐量
ax4b = subplot(2, 2, 2); hold(ax4b, 'on'); grid(ax4b, 'on'); box(ax4b, 'on');
plot(ax4b, snrRange, epThEst, '-o', 'Color', cEpE, 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'MarkerFaceColor', cEpE);
plot(ax4b, snrRange, gfThEst, '-^', 'Color', cGfE, 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'MarkerFaceColor', cGfE);
xlabel(ax4b, 'SNR (dB)', 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel(ax4b, 'bps/Hz', 'FontSize', 11, 'FontName', 'Times New Roman');
title(ax4b, '(b) Throughput (Estimated CSI)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax4b, {'EP','GI-Free'}, 'Location', 'southeast', 'FontSize', 9);
set(ax4b, 'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.5);

% (c) 资源构成
ax4c = subplot(2, 2, 3); hold(ax4c, 'on'); grid(ax4c, 'on'); box(ax4c, 'on');
b4c = bar(ax4c, bm', 'stacked');
for k = 1:4, b4c(k).FaceColor = ovhC(k, :); b4c(k).EdgeColor = 'none'; end
set(ax4c, 'XTickLabel', {'EP','GI-Free'}, 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel(ax4c, 'Samples', 'FontSize', 11, 'FontName', 'Times New Roman');
title(ax4c, '(c) Resource Breakdown', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax4c, {'CPP','ZP','Pilot','Data'}, 'Location', 'north', 'FontSize', 8);
set(ax4c, 'Color', bgC, 'GridColor', gridC, 'GridAlpha', 0.5);
ylim(ax4c, [0 physicalFrameLength * 1.12]);

% (d) 参数表
ax4d = subplot(2, 2, 4); axis(ax4d, 'off');
tbl = {
    '参数', 'EP-AFDM', 'GI-Free';
    '空口采样', num2str(physicalFrameLength), num2str(physicalFrameLength);
    'CPP', num2str(epCfg.PrefixLength), '0';
    'N_data', num2str(epCfg.NumDataSubcarriers), num2str(gfCfg.NumSubcarriers);
    '导频', sprintf('Single (%d)', epCfg.PilotSequenceLength), 'Single (1)';
    'ZP', sprintf('2×%d=%d', epCfg.ZeroPaddingLength, 2*epCfg.ZeroPaddingLength), '0';
    '有效数据', num2str(epDataSymbols), num2str(gfDataSymbols);
    '开销', sprintf('%.1f%%', epOverhead*100), sprintf('%.1f%%', gfOverhead*100);
    '效率', sprintf('%.3f', epSe), sprintf('%.3f', gfSe);
};
numR = size(tbl, 1); numC = size(tbl, 2);
rH = 0.85 / numR;
cW = [0.26 0.37 0.37];
sX = 0; sY = 0.92;
for r = 1:numR
    yP = sY - (r - 1) * rH; xP = sX;
    for c = 1:numC
        if r == 1, fw = 'bold'; bg = [0.88 0.91 0.95];
        else, fw = 'normal'; if mod(r,2)==0, bg = [0.96 0.96 0.96]; else, bg = [1 1 1]; end; end
        rectangle(ax4d, 'Position', [xP yP-rH cW(c) rH], ...
            'FaceColor', bg, 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.5);
        text(ax4d, xP + cW(c)/2, yP - rH/2, tbl{r,c}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'FontWeight', fw, 'FontName', 'Times New Roman');
        xP = xP + cW(c);
    end
end
title(ax4d, '(d) Parameters', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlim(ax4d, [-0.02 1.02]); ylim(ax4d, [sY - numR*rH - 0.05, sY + 0.05]);

sgtitle(fig4, sprintf('EP vs GI-Free: Single Pilot (Frame = %d)', physicalFrameLength), ...
    'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

exportgraphics(fig4, 'fig4_dashboard.png', 'Resolution', 300);
fprintf('已保存: fig4_dashboard.png\n');

%% ================= 保存 =================
save('simulation_results.mat', ...
    'snrRange', 'numTrials', 'physicalFrameLength', ...
    'epBerEst', 'epBerPerf', 'gfBerEst', 'gfBerPerf', 'gfNmse', ...
    'epThEst', 'epThPerf', 'gfThEst', 'gfThPerf', ...
    'epSe', 'gfSe', 'epOverhead', 'gfOverhead', 'seGain');

fprintf('\n====== 完成 ======\n');
fprintf('EP  有效数据: %d / %d\n', epDataSymbols, physicalFrameLength);
fprintf('GF  有效数据: %d / %d\n', gfDataSymbols, physicalFrameLength);
fprintf('频谱效率提升: +%.1f%%\n', seGain);

%% ================= 信道生成器 =================
function [delays, dopplers, gains] = generateChannel(maxDelay, maxDopIdx, nPaths)
    allDelays = 0:maxDelay;
    delays = sort(allDelays(randperm(length(allDelays), nPaths)))';
    randomPhases = -pi + 2 * pi * rand(nPaths, 1);
    dopplers = maxDopIdx * cos(randomPhases);
    gains = sqrt(1 / (2 * nPaths)) * (randn(nPaths, 1) + 1j * randn(nPaths, 1));
end