"""B8-GC 全局图件 spec。每个图脚本顶部 `from style import *; set_style()`。
集中字号/格式/route 语义,让逐图精修保持最小改动。"""
import matplotlib as mpl
import numpy as np
from matplotlib.ticker import MaxNLocator

# ── 规范化 route 顺序(生物学路径)+ 缩写 ──────────────────────
ROUTE_ORDER_BIO = ["primary", "ovarian_met", "ascites", "peritoneal_met",
                   "peritoneal_lavage", "LN_met", "liver_met",
                   "adjacent_normal", "precancer"]
ROUTE_ABBR = {
    "primary": "primary", "ovarian_met": "ovarian", "ascites": "ascites",
    "peritoneal_met": "peri", "peritoneal_lavage": "lavage",
    "LN_met": "LN", "liver_met": "liver",
    "adjacent_normal": "adj-norm", "precancer": "precancer",
}


def abbr(r):
    return ROUTE_ABBR.get(str(r), str(r))


def order_routes(present):
    """按生物学顺序排列出现的 route,未知的追加在后。"""
    present = list(present)
    out = [r for r in ROUTE_ORDER_BIO if r in present]
    out += [r for r in present if r not in ROUTE_ORDER_BIO]
    return out


def set_style():
    """一次性设全局 rcParams。每个图脚本调一次。"""
    mpl.rcParams.update({
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.titleweight": "bold",
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "legend.frameon": False,
        "legend.handletextpad": 0.4,
        "legend.columnspacing": 0.8,
        "legend.borderaxespad": 0.0,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "lines.linewidth": 1.0,
        "figure.dpi": 110,
        "savefig.dpi": 350,
        "savefig.bbox": "tight",      # 关键:外置 legend 不被裁掉
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,           # Illustrator 可编辑文字
        "ps.fonttype": 42,
        "axes.formatter.useoffset": False,   # 关掉 +1e3 偏移
        "axes.formatter.limits": (-4, 6),    # 推迟科学计数法
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    })


def legend_outside(ax, ncol=1, fontsize=6, loc="upper left",
                   anchor=(1.01, 1.0), **kw):
    """把 legend 放到坐标轴右外侧,不盖数据。"""
    return ax.legend(loc=loc, bbox_to_anchor=anchor, ncol=ncol,
                     fontsize=fontsize, frameon=False, borderaxespad=0.0,
                     handletextpad=0.4, columnspacing=0.8, **kw)


def int_ticks(ax, axis="y"):
    """整数刻度(去掉 3.00/2.75 这种脚本感)。"""
    if axis in ("y", "both"):
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    if axis in ("x", "both"):
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))


def robust_vmax(vals, q=98):
    """分位数 vmax,提升 feature/score 图对比度。"""
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, q)) if v.size else None
