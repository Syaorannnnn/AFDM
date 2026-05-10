"""
ID2P 相关径推广 Python 包。

供 MATLAB pyrunfile() 调用或独立运行 runner.py。

模块:
  channel  — BDLR 相关信道 + Dirichlet 核
  id2p     — ID2P 频谱分解 (Delta, 空间协方差)
  cfar     — M_eff 有效候选数 + 门限修正
  runner   — 实验编排 + 可视化
"""

from .channel import (
    compute_chirp_params,
    sample_path_geometry,
    sample_clustered_geometry,
    bdlr_gains,
    build_path_matrix,
    compute_effective_channel,
)

from .id2p import (
    compute_id2p_spectrum,
    compute_effective_id2p,
    compute_null_hotspot_map,
    compute_spatial_cov,
)

from .cfar import (
    compute_effective_candidates,
    compute_independent_threshold,
    compute_corrected_threshold,
    compute_threshold_ratio,
)
