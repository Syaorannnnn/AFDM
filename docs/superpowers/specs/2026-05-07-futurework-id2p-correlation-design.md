# Future Work: Generalized ID2P Modeling Under Correlated Path Gains for AFDM GI-Free Systems

**Date:** 2026-05-07
**Status:** Design spec (implementation pending)
**Scope:** Next publication (journal/conference)

---

## 1. Problem Statement

In the current thesis, the GI-Free AFDM system models path gains as i.i.d. $\mathcal{CN}(0, 1/P)$. Under this assumption, the DAFT-domain ID2P covariance is a scalar matrix $\sigma_w^2 \mathbf{I}_N$, and the CLIP receiver's CFAR threshold uses a single residual power $P_{\mathrm{res}}$ across all candidate positions.

However:
- 3GPP TR 38.901 CDL channel models inherently contain cluster structure with sub-ray merging that produces *de facto* intra-cluster correlation
- V2X, HST, LEO satellite, and ELAA near-field scenarios exhibit non-WSSUS behavior with genuine inter-path gain correlation
- No published AFDM/OTFS work has derived the ID2P covariance in closed form under non-diagonal $\mathbf{R}$

**Core contribution:** Extend the ID2P model from i.i.d. path gains to correlated path gains under a physically grounded cluster correlation structure (BDLR), derive the exact non-scalar ID2P covariance, analyze its subspace structure (coherent cluster contributions produce natural nulls in DAFT domain), and propose CFAR threshold corrections for colored residual interference.

---

## 2. Generalized Channel Model: BDLR

### 2.1 Why Not Weichselberger

**Ruling:** The Weichselberger model $\mathbf{R} = \mathbf{U}_R \tilde{\mathbf{R}} \mathbf{U}_R^H$ is **inappropriate** for single-antenna AFDM.

**Proof (degeneration theorem):** In single-antenna AFDM, the channel is a $P \times 1$ path-gain vector $\mathbf{h}$. With $N_T = 1$, the Weichselberger tensor-product eigenbasis $\mathbf{U}_T^* \otimes \mathbf{U}_R$ collapses to $\mathbf{U}_R$, the coupling matrix $\tilde{\boldsymbol\Omega}$ becomes $P \times 1$, and $\tilde{\mathbf{R}} \equiv \boldsymbol\Lambda$ is strictly diagonal (by definition of $\mathbf{U}_R$ as the eigenbasis of $\mathbf{R}$). The model reduces to the standard KL expansion with zero structural constraint beyond Hermitian PSD.

**Recommendation:** Use KL/eigendecomposition terminology with a footnote explaining the Weichselberger degeneration, then impose physical structure via BDLR.

### 2.2 BDLR Definition

$$\boxed{\mathbf{R} = \mathrm{blkdiag}(\mathbf{R}_1, \dots, \mathbf{R}_C), \qquad \mathbf{R}_c = P_c \rho_c \mathbf{1}\mathbf{1}^H + P_c(1-\rho_c) \mathrm{diag}(p_{c,1}, \dots, p_{c,P_c})}$$

**Parameter mapping from 3GPP TR 38.901 V18 Table 7.7.1:**

| Symbol | Physical Quantity | Source |
|--------|-------------------|--------|
| $C$ | Number of clusters | $= N$ (CDL cluster count) |
| $P_c$ | Rays per cluster | 20 (standard), 1 (LOS cluster) |
| $P_c$ | Cluster total power | Table $P_n$ (dB) → linear |
| $\rho_c \in [0,1]$ | Intra-cluster correlation coefficient | $\rho_c \approx \exp(-\kappa_1 \cdot c_\mathrm{ASA}^2) \cdot \exp(-\kappa_2 \cdot (c_\tau/\tau_\mathrm{rms})^2)$ |
| $\{p_{c,k}\}$ | Sub-ray power weights | $1/P_c$ (uniform) |

**Degeneration checks:**
- $\rho_c = 0$ for all $c$ → $\mathbf{R} = \mathrm{diag}(P_1, \dots, P_P)$, i.i.d. model
- $C = 1$, $P_1 = P$ → full KL expansion
- Rationale: Saleh-Valenzuela (1987), Spencer cluster angular spread (2000), QuaDRiGa implementation (2014)

---

## 3. Exact ID2P Covariance (Novel Contribution)

### 3.1 General Form

DAFT-domain effective channel: $\mathbf{H}_{\mathrm{eff}} = \sum_{i=1}^{P} h_i \mathbf{\Phi}_i$, where $\mathbf{\Phi}_i$ is a unitary operator determined by $(l_i, \nu_i, c_1, c_2, N)$ with $\mathbf{\Phi}_i \mathbf{\Phi}_i^H = \mathbf{I}_N$ (Zheng 2024, Appendix A).

Under $\mathbf{h} \sim \mathcal{CN}(\mathbf{0}, \mathbf{R})$ and i.i.d. white data $\mathbf{x}_d$ with $\mathbb{E}[\mathbf{x}_d \mathbf{x}_d^H] = \sigma_d^2 \mathbf{I}_N$:

$$\boxed{\mathbb{E}_{\mathbf{x}_d}[\mathbf{C}_{\mathrm{ID2P}}] = \sigma_d^2 \sum_{i,j=1}^{P} R_{ij} \mathbf{\Phi}_i \mathbf{\Phi}_j^H + N_0 \mathbf{I}_N}$$

**Degeneration check:** $R_{ij} = \sigma_{h_i}^2 \delta_{ij}$ → $(\sum \sigma_{h_i}^2) \sigma_d^2 \mathbf{I}_N + N_0 \mathbf{I}_N$, matching Zheng 2024 Lemma 2.

### 3.2 BDLR Simplification

Substituting the BDLR structure:

$$\mathbb{E}_{\mathbf{x}_d}[\mathbf{C}_{\mathrm{ID2P}}] = \underbrace{\sigma_d^2 \sum_{c=1}^{C} P_c \rho_c \mathbf{G}_c \mathbf{G}_c^H}_{\text{coherent cluster contributions (rank-1 each)}} + \underbrace{\sigma_d^2 \left(\sum_c P_c(1-\rho_c)\right) \mathbf{I}_N + N_0 \mathbf{I}_N}_{\text{scalar background (white)}}$$

where $\mathbf{G}_c := \sum_{i \in c} \mathbf{\Phi}_i \mathbf{x}_d$ is the coherent cluster kernel in the DAFT domain.

### 3.3 Key Physical Observations

1. **Nulls appear naturally:** When $\rho_c \to 1$ (narrow angular spread, mmWave), cluster $c$ contributes a rank-1 term $P_c \cdot \mathbf{G}_c \mathbf{G}_c^H$. In the nullspace of $\mathbf{G}_c$, ID2P from that cluster **completely cancels**.
2. **Whiteness breaks:** Diagonal entries $[\mathbb{E} \mathbf{C}_{\mathrm{ID2P}}]_{nn}$ are no longer constant — their variance is governed by $|\mathbf{G}_c[n]|^2$ fluctuation.
3. **Computational efficiency:** The double sum of $P^2$ terms compresses to $C$ outer products ($C \ll P$), each $O(N)$.

### 3.4 Literature Positioning

**No published AFDM/OTFS work** (Bemani 2023, Yin 2022/2024, Zheng 2024, Wang 2024, Raviteja 2019, Mishra 2022, Yuan 2021) has derived this non-diagonal ID2P covariance. All existing works assume $\mathbf{R} = \mathrm{diag}$ and obtain scalar $\sigma_w^2 \mathbf{I}_N$. The double-sum closed form is a new independent contribution.

---

## 4. CLIP Receiver Corrections for Colored Residual

### 4.1 Gaussianity of Residual

Under correlated path gains, the residual $\mathbf{e} = \sum_i h_i \mathbf{\Phi}_i \mathbf{x}_d + \mathbf{z}$ remains **complex Gaussian** (linear combination of Gaussian path gains), but its covariance is **colored**:

$$\mathbf{e} \sim \mathcal{CN}(\mathbf{0}, \mathbf{C}_e), \quad \mathbf{C}_e = \sigma_d^2 \sum_{i,j} R_{ij} \mathbf{\Phi}_i \mathbf{\Phi}_j^H + N_0 \mathbf{I}_N$$

**Scalar approximation validity condition:** $\kappa(\mathbf{C}_e) \leq 1.26$ (~1 dB) → threshold error < 1 dB. When $\kappa \leq 2$ → error < 3 dB.

### 4.2 CFAR Threshold Correction

The candidate matching metric $c_k = |\mathbf{h}_k^H \mathbf{r}|^2 / \|\mathbf{h}_k\|^2$ under $\mathcal{H}_0$:

$$\mathbf{h}_k^H \mathbf{e} \sim \mathcal{CN}(0, \mathbf{h}_k^H \mathbf{C}_e \mathbf{h}_k)$$

Define **position-dependent effective residual power** (Rayleigh quotient):

$$P_{\mathrm{res}}^{\mathrm{eff}}(k) := \frac{\mathbf{h}_k^H \mathbf{C}_e \mathbf{h}_k}{\|\mathbf{h}_k\|^2} \in [\lambda_{\min}(\mathbf{C}_e), \lambda_{\max}(\mathbf{C}_e)]$$

**Corrected CFAR threshold (per-candidate):**

$$\eta_k^{\mathrm{corr}} = s^{(k)} \cdot P_{\mathrm{res}}^{\mathrm{eff}}(k) \cdot \ln\left(\frac{M_{\mathrm{cand}}}{P_{\mathrm{FA}}}\right)$$

**Corrected CFAR threshold (conservative scalar):**

$$\eta_{\mathrm{CFAR}}^{\mathrm{corr}} = s^{(k)} \cdot P_{\mathrm{res}}^{\max} \cdot \ln\left(\frac{M_{\mathrm{cand}}}{P_{\mathrm{FA}}}\right), \quad P_{\mathrm{res}}^{\max} = \max_k P_{\mathrm{res}}^{\mathrm{eff}}(k)$$

**Degeneration check:** When $\mathbf{R}$ is diagonal, $P_{\mathrm{res}}^{\mathrm{eff}}(k) = P_{\mathrm{res}}$ for all $k$, recovering the original formula.

### 4.3 Progressive Strictness Factor Extension

Three correction levels (increasing complexity):

1. **Conservative scalar gain:** $s_{\mathrm{corr}}^{(k)} = s^{(k)} \cdot \sqrt{\lambda_{\max}(\tilde{\mathbf{R}}) / \lambda_{\mathrm{avg}}(\tilde{\mathbf{R}})}$, where $\tilde{\mathbf{R}} = \mathbf{R} / (\mathrm{Tr}(\mathbf{R})/P)$
2. **Per-candidate adaptive:** $\eta_k = s^{(k)} \cdot P_{\mathrm{res}}^{\mathrm{eff}}(k) \cdot \ln(M_{\mathrm{cand}}/P_{\mathrm{FA}})$ (aligned with Robey AMF)
3. **Condition-number margin:** When $\kappa(\mathbf{C}_e) > \kappa_{\mathrm{th}}$, add $\zeta \cdot \log(\kappa(\mathbf{C}_e))$, $\zeta \in [0.1, 0.3]$

### 4.4 Safety Property (Non-Intuitive)

Positive correlation among candidate metrics $c_k$ makes the original Bonferroni union bound **conservative** (Das-Bhandari 2019). This means: **the existing CFAR threshold under correlated paths does not increase false alarm probability — it only loses detection probability.** The corrections above recover lost $P_d$ without risk of increasing $P_{\mathrm{FA}}$.

---

## 5. Simulation Verification Plan

### 5.1 Channel Configuration

Generate channels using QuaDRiGa or NYUSIM with 3GPP CDL-C/D parameters:

| Parameter | Value |
|-----------|-------|
| Carrier frequency | 4 GHz |
| Bandwidth | 5 MHz ($N = 512$, $\Delta f \approx 1$ kHz) |
| CDL profile | C (NLOS, 24 clusters) or D (LOS, 13 clusters) |
| Per-cluster rays | 20 (sub-ray merging to cluster taps) |
| $\rho_c$ sweep | {0, 0.3, 0.6, 0.9} |
| SNR range | 0–25 dB |

### 5.2 Baseline Comparisons

| Configuration | Description |
|---------------|-------------|
| A | i.i.d. $\mathbf{h}$ + original CLIP (thesis baseline) |
| B | BDLR $\mathbf{h}$ + original CLIP (degradation measurement) |
| C | BDLR $\mathbf{h}$ + scalar-corrected CFAR (level 1) |
| D | BDLR $\mathbf{h}$ + per-candidate adaptive CFAR (level 2) |

### 5.3 Metrics

- BER, $P_d$, NMSE (standard)
- $\kappa(\hat{\mathbf{C}}_e)$ distribution over iterations
- Per-position $P_{\mathrm{res}}^{\mathrm{eff}}(k)$ histogram
- Null depth: $\min_n |\mathbf{G}_c[n]|^2$ vs $\rho_c$

### 5.4 Validation Steps

1. **Degeneration check:** Verify $\mathbf{C}_{\mathrm{ID2P}}$ → scalar when $\rho_c = 0$
2. **Empirical $\hat{\mathbf{R}}$ estimation:** From QuaDRiGa realizations, compute Frobenius error $\|\hat{\mathbf{R}} - \hat{\mathbf{R}}_{\mathrm{BDLR}}\|_F$
3. **Bures distance:** Between empirical and BDLR-fitted covariance for ID2P
4. **CFAR validity:** Monte Carlo empirical $P_{\mathrm{FA}}$ vs target under BDLR with each correction level

---

## 6. Thesis Chapter 3 Integration Points

Five revision locations identified (referencing thesis structure):

1. **§3.1 (Channel model):** Add paragraph acknowledging 3GPP Uncorrelated Scattering default, then generalize to $\mathbf{h} \sim \mathcal{CN}(\mathbf{0}, \mathbf{R})$
2. **§3.2 (ID2P derivation):** Add section "Exact Covariance Under Path Correlation" ($\S$4.2 of this doc), with degeneration-check footnote
3. **§3.3 (BDLR):** Introduce cluster block-diagonal structure with 3GPP parameter mapping — use KL terminology, not Weichselberger
4. **§3.4 (CLIP receiver):** Add "Colored Residual Correction" subsection covering modified CFAR threshold and progressive strictness extension
5. **§3.5 (Complexity note):** $O(NC)$ vs $O(NP^2)$ for BDLR

---

## 7. Key References

**3GPP & Channel Modeling:**
- 3GPP TR 38.901 V19.0.0 (2025-06)
- Saleh & Valenzuela, IEEE JSAC 1987
- Spencer et al., IEEE JSAC 2000
- Jaeckel et al. (QuaDRiGa), IEEE TAP 2014
- Bernadó et al., IEEE TVT 2014 (V2V non-stationary)

**Weichselberger & MIMO Correlation:**
- Weichselberger et al., IEEE TWC 2006
- Sayeed, IEEE TSP 2002
- Özcelik et al., Electron. Lett. 2003 (Kronecker deficiency)

**AFDM / OTFS:**
- Bemani et al., IEEE TWC 2023 (AFDM foundation)
- Zheng et al., arXiv:2404.10232 (superimposed pilots)
- Yin et al., IEEE/CIC ICCC 2022; IEEE TWC 2024
- Wang et al., arXiv:2404.01088 (GI-Free)
- Raviteja et al., IEEE TVT 2019 (OTFS EP)
- Mishra et al., IEEE TWC 2022 (OTFS superimposed)

**CFAR & Colored Detection:**
- Kelly, IEEE TAES 1986 (GLRT)
- Robey et al., IEEE TAES 1992 (AMF)
- Das & Bhandari, arXiv:1908.02193 (FWER under correlation)
- Nadakuditi & Silverstein, IEEE JSTSP 2010 (condition number)

---

## 8. Open Questions

1. **BDLR parameter estimation from real channels:** Can $\rho_c$ be reliably estimated from a limited number of QuaDRiGa realizations? Minimum sample complexity?
2. **$\kappa(\mathbf{C}_e)$ typical values:** For CDL-C/D at 4 GHz, what is the empirical distribution of $\kappa(\mathbf{C}_e)$ across SNR? Does it cross the 1.26 threshold in practical regimes?
3. **Null depth exploitation:** Can the natural DAFT-domain nulls from coherent clusters be exploited for *active* data precoding (placing data symbols in null directions), rather than just passive CFAR correction?
4. **Extension to MIMO-AFDM:** Does the BDLR structure naturally extend to multi-antenna via $\mathbf{R}_{\mathrm{BS}} \otimes \mathbf{R}_{\mathrm{cluster}}$ Kronecker form, or is a more general jointly-correlated structure required?
