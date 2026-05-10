"""
matplotlib 全局配置：CJK 字体支持 + 样式预置。

所有 py_correlated 绘图脚本应在 import matplotlib.pyplot 前
调用 setup_plots()。
"""

import matplotlib
import matplotlib.pyplot as plt


def setup_plots():
    """配置 matplotlib 使用中文字体，禁用 Type 3 字体。"""
    matplotlib.rcParams.update({
        'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'SimSun',
                            'DejaVu Sans'],
        'axes.unicode_minus': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })
    # 清除字体缓存确保新字体生效
    matplotlib.font_manager._load_fontmanager(try_read_cache=False)
