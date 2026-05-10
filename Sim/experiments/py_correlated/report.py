"""
Python 级验证报告生成器。

输出 Markdown 摘要，记录配置、ID2P 数据子载波维度 M_eff、
候选网格维度 M_eff、经验 Pfa/Pd 和主要解释边界。
"""

from pathlib import Path


def write_markdown_report(path, cfg, id2p_rows, grid_rows):
    """
    写入 Markdown 验证报告。

    参数:
      path:       输出文件路径
      cfg:        实验配置字典
      id2p_rows:  list[dict] ID2P 子载波维度 rho 扫描行，
                  每行含 rho / m_eff / threshold_ratio
      grid_rows:  list[dict] 候选网格维度验证行，
                  每行含 rho / m_eff_grid / empirical_pfa_* / empirical_pd_*
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Python CFAR M_eff Validation Report",
        "",
        "## Configuration",
        "",
    ]
    for key in sorted(cfg):
        lines.append(f"- `{key}`: `{cfg[key]}`")

    # ID2P 子载波维度
    lines.extend(["", "## ID2P Subcarrier Dimension", ""])
    if id2p_rows:
        lines.append("| rho | M_eff | threshold_ratio |")
        lines.append("|---:|---:|---:|")
        for row in id2p_rows:
            lines.append(
                f"| {row['rho']:.2f} | {row['m_eff']:.3f} | "
                f"{row['threshold_ratio']:.4f} |"
            )
    else:
        lines.append("(no data)")

    # 候选网格维度
    lines.extend(["", "## Candidate Grid Dimension", ""])
    if grid_rows:
        lines.append(
            "| rho | dir_r | M_eff_grid | M_raw | ratio | "
            "Pfa_ind | Pfa_corr | Pd_ind | Pd_corr |"
        )
        lines.append(
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
        )
        for row in grid_rows:
            dr = row.get("dirichlet_r", "")
            lines.append(
                f"| {row.get('rho', '?'):.2f} | {dr} | "
                f"{row.get('m_eff_grid', 0):.3f} | "
                f"{row.get('m_raw_grid', 0):.0f} | "
                f"{row.get('threshold_ratio', 0):.4f} | "
                f"{row.get('empirical_pfa_independent', 0):.4g} | "
                f"{row.get('empirical_pfa_corrected', 0):.4g} | "
                f"{row.get('empirical_pd_independent', 0):.4g} | "
                f"{row.get('empirical_pd_corrected', 0):.4g} |"
            )
    else:
        lines.append("(no data)")

    # 解释规则
    lines.extend([
        "",
        "## Interpretation Boundary",
        "",
        "ID2P 子载波维度的 M_eff 用于解释数据泄漏场相关性。",
        "候选网格维度的 M_eff 才能作为后续 CFAR 门限公式调整的直接依据。",
        "本报告不调用 MATLAB 主干仿真，因此不直接声明 BER 改善。",
        "",
        "## Decision Rules",
        "",
        "- 若 ID2P 子载波维度 M_eff 很小，但候选网格维度 M_eff 不小，"
        "则只说明数据泄漏场相关，不支持直接调整 OMP-CFAR 候选数。",
        "- 若候选网格维度 M_eff 也显著小于原始候选数，且修正门限提高 Pd，"
        "同时经验 Pfa 未失控，则该修正可进入 MATLAB 端候选方案。",
        "- 若修正门限提高 Pd 但经验 Pfa 明显高于目标 Pfa，"
        "应保留为诊断结论，不进入主接收机。",
        "- 本阶段不声明 BER 改善，因为 Python 验证没有执行完整数据检测链路。",
    ])

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(out_path)
