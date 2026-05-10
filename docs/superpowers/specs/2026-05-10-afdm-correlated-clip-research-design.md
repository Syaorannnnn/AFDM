# AFDM-GI-Free 相关径 CLIP 学硕研究方向：理论+算法 Design Spec

**Date:** 2026-05-10
**Status:** Research design (brainstorming approved, not yet implementation)
**Scope:** 学硕研究周期内的 primary research direction，目标 TWC/TSP 级 journal + 前置 letter
**Supersedes / Extends:** `2026-05-07-futurework-id2p-correlation-design.md`（该文档覆盖本 spec 中 M1-B §3.1 部分，本 spec 扩展其余 M1-B、M1-C、M2、评估、时间线）

---

## 1. 研究定位

**Primary axis (M1, 理论)**：首次为 AFDM-GI-Free 在相关径（BDLR）信道下推导 ID2P 的严格统计模型，并据此给出 **LRT-最优路径检测器、参数估计 CRLB、与新 CFAR 门限**。理论深度停在 "B 层"——四个闭式 + 一个条件数 $\kappa(\boldsymbol{\Sigma})$ 上界，不进 random matrix 渐近。

**Secondary axis (M2, 算法)**：在 M1-B 派生的 Capon-style generalized matched filter 骨架上，用 **unrolled EM network (M2-α)** 学习 $\boldsymbol{\Sigma}$ 超参，辅以 **门限校准模块 (γ)**，端到端优化 path detection + channel estimation。可选 **contrastive cluster assignment (β)** 作为 extension。

**Stretch goal (M1-C)**：**可辨识性 phase transition** 的数值 evidence + 闭式猜想，作为 appendix 或 future work。不押注成败。

**Validation environment (M3, 辅)**：在 3GPP CDL-C/D 信道下提供系统性 benchmark（vs OTFS-CLIP、OFDM-DMRS、AFDM-EP）。M3 作为实验章节的 validation field，不独立成 axis。

**论文叙事链**：correlated multipath → identifiability crisis → LRT-CRLB 理论根基 → unrolled ML 算法 → 四件套 metric 验证 → 3GPP CDL benchmark。

**比例**：按 journal 25 页估计，通信原理约 50%，机器学习约 30%，评估与验证约 20%。审稿买单的 **key translator** 是 $\kappa(\boldsymbol{\Sigma})$ 上界——它既是 M1-B 理论收尾，又是 M2 的 L3 loss 数值约束，又是 M3 metric 的归一化基线，三处引用让"通信+ML"真正合为一篇论文。

---

## 2. 与 2026-05-07 先行设计的关系

`2026-05-07-futurework-id2p-correlation-design.md` 已完成：

- BDLR 信道模型定义（§2.2）与 Weichselberger degeneration 论证（§2.1）
- ID2P 协方差双 sum 闭式（§3.1）及 BDLR 简化（§3.2）
- Position-dependent CFAR 门限 $P_{\mathrm{res}}^{\mathrm{eff}}(k)$ 与保守标量版本（§4.2）
- Progressive strictness 扩展三层（§4.3）
- Chapter 3 integration plan（§6）

本 spec **不重复**以上内容，直接 cite。本 spec 在此基础上新增：**LRT detector、CRLB、$\kappa$ 上界（M1-B 扩展）**、**M1-C phase transition**、**M2-α 算法与训练**、**四件套评估**、**学硕时间线**、**Python 主线工程规划**。

---

## 3. Section M1-B 扩展：LRT Detector、CRLB 与 $\kappa$ 上界

### 3.1 LRT-最优 path detector（novel contribution）

2026-05-07 §4.2 给出的 $P_{\mathrm{res}}^{\mathrm{eff}}(k)$ 是对标量 CFAR 的位置相关修正；本节进一步给出**LRT-最优 detector 的完整结构**。

记 $\boldsymbol{\Sigma} = \mathbf{C}_{\mathrm{ID2P}} + N_0 \mathbf{I}_N$（位置-色化残差协方差，含未知路径 ID2P + 噪声）。对候选 $k$：

$$\mathcal{H}_0: \mathbf{r} = \mathbf{i}_{\mathrm{other}} + \mathbf{z}, \quad \mathcal{H}_1: \mathbf{r} = h_k \boldsymbol{\Phi}_k \mathbf{x}_p + \mathbf{i}_{\mathrm{other}} + \mathbf{z}$$

Generalized Likelihood Ratio Test 在 $h_k$ 未知 + $\boldsymbol{\Sigma}$ 已知下导出 **Capon-style matched filter**：

$$\boxed{T(\mathbf{r}; k) = \frac{|\boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma}^{-1} \mathbf{r}|^2}{\boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k}}$$

对比现有 CLIP 代码 `group_omp.py::group_omp_detect` 中的 naive matched filter $|\boldsymbol{\Phi}_k^{\mathsf H} \mathbf{r}|^2$，差别在于 $\boldsymbol{\Sigma}^{-1}$ 白化。SNR 增益定义为

$$\mathrm{SNR}_{\mathrm{gain}}(k) = 10 \log_{10}\left(\frac{\boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k \cdot \boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma} \boldsymbol{\Phi}_k}{\|\boldsymbol{\Phi}_k\|^4}\right)$$

当 $\boldsymbol{\Sigma} \propto \mathbf{I}_N$（i.i.d. 极限）该增益为 0 dB，验证 degeneration。当 $\rho \to 1$ 时 $\boldsymbol{\Sigma}$ rank 下降，$\boldsymbol{\Sigma}^{-1}$ 需用 pseudo-inverse 或 Tikhonov，此时增益可达 5–10 dB，是本论文的主要性能来源。

### 3.2 Fisher 信息矩阵与 CRLB

对参数 $\boldsymbol{\theta} = [\tau, \alpha, \mathrm{Re}\,h, \mathrm{Im}\,h]^{\mathsf T}$（单径 delay、Doppler、增益实部虚部）：

$$J_{ij}(\boldsymbol{\theta}) = 2\,\mathrm{Re}\!\left[\!\left(\frac{\partial\boldsymbol{\mu}}{\partial\theta_i}\right)^{\!\mathsf H}\!\boldsymbol{\Sigma}^{-1}\!\left(\frac{\partial\boldsymbol{\mu}}{\partial\theta_j}\right)\!\right], \qquad \boldsymbol{\mu}(\boldsymbol{\theta}) = h \boldsymbol{\Phi}(\tau,\alpha) \mathbf{x}_p$$

**关键定理（待推）**：簇内两条路径 $i, j \in c$ 满足 $\cos(\angle(\boldsymbol{\Phi}_i, \boldsymbol{\Phi}_j)) \to 1$ 时（工程实测 0.998），$\mathbf{J}$ 关于这两径的 block 接近奇异，CRLB $[\mathbf{J}^{-1}]_{ii} \to \infty$，**严格证明单径局部不可辨**。

**簇级 CRLB**：定义簇级参数 $\boldsymbol{\theta}_c = [\bar{\tau}_c, \bar{\alpha}_c, \mathrm{Re}\,g_c, \mathrm{Im}\,g_c]^{\mathsf T}$（簇质心 delay、Doppler、聚合增益），簇级 Fisher 矩阵 $\mathbf{J}_{\mathrm{cluster}}$ 对所有 $\rho \in [0,1]$ 满秩，对应 CRLB 有限。这一对比量化了 M1-B 的核心声明"簇化恢复可辨识性"。

### 3.3 $\boldsymbol{\Sigma}$ 条件数上界（Theory→ML translator）

由 2026-05-07 §3.2 的 BDLR 简化 ID2P 协方差

$$\boldsymbol{\Sigma} = \sigma_d^2 \sum_{c=1}^{C} P_c \rho_c \mathbf{G}_c \mathbf{G}_c^{\mathsf H} + \sigma_d^2 \Big(\sum_c P_c(1-\rho_c)\Big) \mathbf{I}_N + N_0 \mathbf{I}_N$$

上界推导（via Weyl 不等式）：

$$\kappa(\boldsymbol{\Sigma}) \leq \frac{\sigma_d^2 \max_c P_c \rho_c \|\mathbf{G}_c\|^2 + \sigma_d^2 \bar{P}(1-\bar{\rho}) + N_0}{\sigma_d^2 \bar{P}(1-\bar{\rho}) + N_0}$$

其中 $\bar{P} = \frac{1}{C}\sum P_c$，$\bar{\rho} = \frac{1}{C}\sum \rho_c$。当 $\rho \to 1$ 时分母 $\to N_0$，条件数以 $\mathrm{SNR}_{\mathrm{p}}$ 速率发散。

**translator 作用**：

- 给 M1-B §4.1（scalar approximation 有效性）提供 analytical 而非 empirical 的 $\kappa$ 阈值。
- 给 M2-α 的 L3 ELBO loss 提供数值稳定性约束：规定训练时 $\hat{\boldsymbol{\Sigma}}$ 必须满足 $\kappa \leq \kappa^\ast$（由论文 Theorem 给出）才计算 $\log\det$，否则 fallback 到 Tikhonov regularized 版本。
- 给 M3 的 Effective DoF metric 提供归一化基线：$\eta_{\mathrm{DoF}}(\boldsymbol{\Sigma}) := \mathrm{rank}_\epsilon(\boldsymbol{\Sigma}; \epsilon = \lambda_{\max}/\kappa^\ast)$。

### 3.4 M1-B Deliverables Checklist

- [ ] D1: ID2P 二阶统计量闭式（已在 2026-05-07 §3 完成）
- [ ] D2: LRT-optimal $T(\mathbf{r}; k)$ 推导 + SNR gain 分析
- [ ] D3: Fisher 矩阵 + 单径 CRLB + 簇级 CRLB 对比定理
- [ ] D4: 新 CFAR 门限 $\eta_{\mathrm{new}} = \ln(M/P_{\mathrm{FA}}) \cdot \boldsymbol{\Phi}_k^{\mathsf H} \boldsymbol{\Sigma}^{-1} \boldsymbol{\Phi}_k \cdot \sigma_{\mathrm{eff}}^2$
- [ ] D5: $\kappa(\boldsymbol{\Sigma})$ 上界与 $\kappa^\ast$ 阈值（translator）
- [ ] D6: Letter 稿（SPL/WCL 级别，基于 D1-D5）

---

## 4. Section M1-C（Stretch Goal）：可辨识性 Phase Transition

**目标**：证明或给出数值 evidence，存在临界曲线 $\rho_c^\ast(\mathrm{SNR}, P_c, \alpha_{\max})$，使

- $\rho > \rho_c^\ast$：$\mathbf{J}$ 的 $\epsilon$-rank 塌缩，单径局部不可辨；
- $\rho < \rho_c^\ast$：$\mathbf{J}$ 满秩，单径局部可辨；
- $\forall \rho$：簇级 $\mathbf{J}_{\mathrm{cluster}}$ 满秩（簇化恢复定理）。

**工具路线**（优先级从低到高风险）：

1. **Concentration + large deviation**：给出 $\mathbb{P}[\mathrm{rank}(\mathbf{J}) < P]$ 的指数衰减界，临界 $\rho$ 由衰减速率等于 0 处定义。
2. **Random matrix 渐近**（Marchenko-Pastur 风格）：$P \to \infty, P_c = \Theta(P^\beta)$ scaling 下奇异值 bulk 与 0 碰撞点。
3. **Information-theoretic**：互信息 $I(\boldsymbol{\theta}; \mathbf{y})$ 的 saturation 刻画。

**降级策略**：若推导不出闭式 $\rho_c^\ast$，降为 numerical observation 章节，给出 3GPP CDL-C/D 下 $\hat{\rho}_c^\ast$ 的 Monte Carlo 拟合曲线 + 开放猜想。不影响 M1-B 的完整性。

---

## 5. Section M2-α：Unrolled EM Network + γ 门限校准

### 5.1 网络结构

Unrolled EM network，共 $K$ 层（初拟 $K = 10$）。每层对应一次完整 EM 迭代：

**E-step (unrolled)**：
$$\hat{\boldsymbol{\Sigma}}_k^{-1} \leftarrow \hat{\boldsymbol{\Sigma}}_{k-1}^{-1} + \eta_k \cdot \mathrm{Neumann}_L(\hat{\boldsymbol{\Sigma}}_{k-1}, \Delta \hat{\boldsymbol{\Sigma}}_k)$$

Neumann 展开阶数 $L = 2$（$\mathbf{A}^{-1}(\mathbf{I} + \Delta\mathbf{A}\mathbf{A}^{-1})^{-1} \approx \mathbf{A}^{-1}(\mathbf{I} - \Delta\mathbf{A}\mathbf{A}^{-1})$），step size $\eta_k$ learnable。

**M-step (differentiable moment matching)**：
$$(\hat{\rho}_c, \hat{\sigma}_c^2)_k \leftarrow \mathrm{MLP}_\theta\!\left(\mathrm{vec}(\hat{\boldsymbol{\Sigma}}_k), \hat{\mathcal{S}}_{k-1}\right)$$

输出受 BDLR 结构约束（simplex / positive orthant projection）。

**LRT 层**：
$$\hat{T}_k(\boldsymbol{\Phi}_j) = \frac{|\boldsymbol{\Phi}_j^{\mathsf H} \hat{\boldsymbol{\Sigma}}_k^{-1} \mathbf{r}|^2}{\boldsymbol{\Phi}_j^{\mathsf H} \hat{\boldsymbol{\Sigma}}_k^{-1} \boldsymbol{\Phi}_j}$$

通过 Sinkhorn 软 top-$K$ 得到软支撑 $\hat{\mathbf{p}}_k \in [0,1]^{M_{\mathrm{cand}}}$。

### 5.2 γ 门限校准模块

最后一层后接 2-layer MLP，输入 $(\hat{\kappa}(\boldsymbol{\Sigma}_K), \log P_{\mathrm{FA}}, \mathrm{SNR})$，输出 scale $s_\gamma \in [0.5, 2.0]$：

$$\tilde{\eta} = s_\gamma(\hat{\kappa}, P_{\mathrm{FA}}, \mathrm{SNR}) \cdot \eta_{\mathrm{new}}$$

训练时约束 $s_\gamma$ 在 $\hat{\kappa} \to 1$ 处收敛到 1（degeneration guarantee）。

### 5.3 β extension（可选）

Contrastive 簇嵌入 + differentiable k-means 学簇分配。只在 T4 阶段有余裕时加入，不进核心 deliverables。

### 5.4 训练方案（D3 curriculum + L5 composite）

**Stage A (月 7-8, D1 全监督)**：
$$\mathcal{L} = 0.7 \mathcal{L}_{\mathrm{supp}} + 0.3 \mathcal{L}_{\mathrm{nmse}}$$

3GPP CDL-C/D + QuaDRiGa 生成 $10^5$ 样本，ground-truth $\{\mathbf{h}_i, \boldsymbol{\Sigma}^\ast, \mathcal{S}^\ast\}$ 全部可得。

**Stage B (月 8-9, D1 加下游对齐)**：
$$\mathcal{L} = 0.5 \mathcal{L}_{\mathrm{supp}} + 0.3 \mathcal{L}_{\mathrm{nmse}} + 0.2 \mathcal{L}_{\mathrm{ber}}$$

从 Stage A checkpoint continue，$\lambda_{\mathrm{ber}}$ 按 epoch linearly 从 0 升到 0.3，避免 BER gradient 震荡。

**Stage C (月 9-10, D2 伪标签 self-supervised)**：
$$\mathcal{L} = 0.5 \mathcal{L}_{\mathrm{ml}} + 0.5 \mathcal{L}_{\mathrm{ber}}$$

teacher = M1-B 的 closed-form LRT detector（无训练），student = unrolled network。L3 ELBO 的 $\log\det$ 项按 §3.3 的 $\kappa^\ast$ 约束 clip。

**Ablation 轴**（自然形成，最大化审稿买单度）：

- A-only / A+B / A+B+C：展示 curriculum 必要性。
- L5 vs L1+L2 / L3-only：展示 composite loss 必要性。
- Unrolled depth $K = \{3, 5, 10, 15\}$：展示理论-网络对齐度。

### 5.5 M2-α Deliverables Checklist

- [ ] E1: PyTorch 实现 unrolled EM K-层网络
- [ ] E2: QuaDRiGa + 3GPP CDL 合成训练集（D1）
- [ ] E3: teacher-student 伪标签管道（D2）
- [ ] E4: L5 composite loss + curriculum scheduler
- [ ] E5: γ 门限校准 MLP
- [ ] E6: Ablation 完整表

---

## 6. 四件套评估指标

保留 **BER** 作为主指标，新增三项 AFDM-GI-Free-correlated 场景**命名化**指标：

1. **Cluster-level $P_d$**：簇支撑命中率 $P_d^{(c)} = |\hat{\mathcal{S}}_c \cap \mathcal{S}_c^\ast| / |\mathcal{S}_c^\ast|$。
2. **Effective DoF ratio**：$\eta_{\mathrm{DoF}} = \mathrm{rank}_\epsilon(\hat{\boldsymbol{\Sigma}}) / \mathrm{rank}_\epsilon(\boldsymbol{\Sigma}^\ast)$，$\epsilon = \lambda_{\max}/\kappa^\ast$。直连 M1-B 定理。
3. **Support Hausdorff distance**：$d_H(\hat{\mathcal{S}}, \mathcal{S}^\ast)$ in $(l, \alpha)$-grid with fractional Doppler scaling.
4. **Column-space alignment**：$\mathrm{align}(\hat{\mathbf{H}}, \mathbf{H}^\ast) = \|\hat{\mathbf{U}}_r^{\mathsf H} \mathbf{U}_r^\ast\|_F^2 / r$，$r = \mathrm{rank}_\epsilon(\mathbf{H}^\ast)$。

**对比矩阵**：

| 轴 | 取值 |
| --- | --- |
| $\rho$ | 0.0 / 0.3 / 0.5 / 0.7 / 0.9 |
| SNR (dB) | 0, 5, 10, 15, 20, 25 |
| 信道 | Jakes i.i.d. / BDLR 合成 / 3GPP CDL-A / CDL-C / CDL-D |
| 方法 | OFDM-DMRS / AFDM-EP / AFDM-GI-Free CLIP (baseline) / LRT-CLIP (M1-B only) / Unrolled-LRT-CLIP (M1-B + M2-α) |

Main text 放 3 个 representative 切片，其余进 supplementary。

---

## 7. 时间线（24 月学硕周期）

| 阶段 | 月 | 核心产出 | 风险缓冲 |
| --- | --- | --- | --- |
| T0 理论准备 | 1–2 | 精读 1-4 + Capon/LRT/CRLB 复习 + 阶 A 工程代码实际度审计 | 1 月 buffer |
| T1 M1-B 推导 | 3–6 | D1-D5 全部完成，D6 letter 投稿 | 若卡 CRLB 块 1 月 |
| T2 M2-α 原型 | 7–10 | E1-E6 全部完成，BER 超越 baseline | 若训练不稳 curriculum 调整 |
| T3 评估建立 | 11–13 | 四件套指标实现 + 3GPP CDL 接口 | QuaDRiGa 接入 2 周 buffer |
| T4 扩展实验 | 14–18 | 完整对比矩阵 + 可选 β extension + M1-C 数值 evidence | β 可砍 |
| T5 论文写作 | 19–22 | Journal 主稿投稿 | 1 月 |
| T6 答辩缓冲 | 23–24 | 修订 + 录用信 + 答辩材料 | |

**Go/no-go gates**：

- T1 底：若 D1-D5 未全齐，letter 推迟，journal 时间线整体延后 1 月。
- T2 底：若 Stage A 训练不收敛，降级到 M2-γ-only（仅门限校准），M2-α 成 appendix。
- T4 中：若 3GPP CDL 接入卡住，降级到 BDLR + Jakes 混合 validation。

---

## 8. Python 主干工程规划

**主路径**：`Sim/experiments/py_correlated/` 继续作为 research codebase。

**新增模块**（估算）：

| 模块 | 路径 | 行数 | 依赖 |
| --- | --- | --- | --- |
| LRT detector (closed-form) | `clip_lrt.py` | ~200 | `grid_stats.py`, `cfar.py` |
| $\boldsymbol{\Sigma}$ 估计与 $\kappa$ 诊断 | `sigma_estimator.py` | ~150 | NumPy + SciPy |
| Unrolled EM network | `unrolled_em.py` | ~400 | PyTorch |
| D1 QuaDRiGa 合成 | `dataset_cdl.py` | ~300 | QuaDRiGa MATLAB bridge or NYUSIM Python |
| D2 teacher-student 管道 | `pseudo_label.py` | ~250 | |
| L5 composite loss | `losses.py` | ~200 | |
| 四件套 metric | `metrics_correlated.py` | ~300 | |
| 对比实验脚本 | `exp_correlated_benchmark.py` | ~400 | |

**MATLAB 主干**：保持当前本科论文 baseline 不动。`matlab_bridge.py` 单向输出 $\{\hat{\rho}, \hat{\sigma}_c^2, \hat{\mathcal{S}}, \eta_{\mathrm{DoF}}\}$ 用于对比图。

**测试覆盖**：每个模块配 pytest unit test；关键 deliverable (D1-D5, E1-E6) 配 golden file 回归测试。

---

## 9. 风险清单

| 风险 | 严重度 | 对策 |
| --- | :---: | --- |
| M1-B CRLB 推导陷入 random matrix | 中 | 停在 finite-$P$ closed-form，必要时引 Gershgorin 界粗略 bound |
| $\boldsymbol{\Sigma}$ 在 $\rho \to 1$ 近奇异 | 中 | §3.3 的 $\kappa^\ast$ 约束显式给定，论文单独一节讨论 |
| M2-α unrolled 训练不稳 | 中 | curriculum + gradient clipping + Stage A checkpoint 保底 |
| D3 伪标签 drift | 低 | Stage C 监控 $\|\hat{\boldsymbol{\Sigma}} - \hat{\boldsymbol{\Sigma}}_{\mathrm{teacher}}\|_F$，超阈回退 |
| QuaDRiGa 接入超预期 | 中 | T3 预留 2 周，卡住切 NYUSIM Python 或 fallback BDLR 合成 |
| M1-C 无闭式 | 低 | 降级 numerical evidence + conjecture，不影响 B |
| T4 被压缩 | 中 | β extension 砍掉，核心 deliverable T3 末须齐 |
| 硕士导师要求回 MATLAB | 低 | 只在论文接受后做 port，不阻塞研究主线 |

---

## 10. 开放问题（from 2026-05-07 §8 + 新增）

承 2026-05-07 §8：

1. BDLR 参数 $\rho_c$ 从有限 QuaDRiGa realizations 的可靠估计样本复杂度。
2. 3GPP CDL-C/D 下 $\kappa(\boldsymbol{\Sigma})$ 的经验分布，是否跨 $\kappa^\ast$ 阈值。
3. Null depth 的主动利用（数据 precoding placing in null space）。
4. MIMO-AFDM 扩展下 BDLR 是否走 Kronecker $\mathbf{R}_{\mathrm{BS}} \otimes \mathbf{R}_{\mathrm{cluster}}$。

新增：

1. $\kappa^\ast$ 阈值的最优选择（CRLB 敏感度 vs ELBO 数值稳定性 trade-off）。
2. Unrolled depth $K$ 与 M1-B 理论 EM 迭代数的 correspondence（是否存在 $K^\ast(\rho, \mathrm{SNR})$ 闭式）。
3. 簇学习（β extension）是否能从 attention 机制直接 induce，而非显式 contrastive。
4. M1-C phase transition 是否可通过 PAC-Bayes 给出 finite-sample 版本。

---

## 11. Deliverable 汇总

**Theoretical (M1-B + M1-C)**：D1-D5 closed-form + D6 letter + M1-C numerical evidence。
**Algorithmic (M2-α + γ)**：E1-E6 PyTorch 实现 + ablation 完整表。
**Metric (四件套)**：M1-M4 命名化指标实现 + 对比矩阵图 + supplementary。
**Paper**：Letter (T1 底) + Journal main (T5) + 可选 M1-C appendix。
**Code**：Python 主干 open source（MIT 或 BSD-3）。

---

## 12. 与本科论文 / memory 的对齐

本 spec 的 M1-B / M1-C / M2-α 均属本科论文 Chapter 5 (Future Work) 第 1 项（"Complexity reduction + theoretical extension"）的 expansion。不修改 `docs/memory/*`，memory 记录本科论文完成的事实，不混入学硕研究设计。

**后续 memory 更新建议**：若 T1 letter 投稿成功，新增 `docs/memory/2026-xx-letter-lrt-clip.md` 记录投稿事实；若 T5 journal 接受，新增 `docs/memory/2026-xx-journal-unrolled-clip.md`。这两项是 memory 的"成果记录"用途，与研究设计文档分离。

