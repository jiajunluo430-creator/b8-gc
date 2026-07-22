"""作图统一样式 + 配色 + 保存。npj Precision Oncology / Nature Portfolio 规范。
Journal of Translational Medicine (BMC): 多panel合成单图, 小写 a/b/c,
- 宽 单栏85mm / 双栏170mm; 最高 ~225mm; 竖版优先(2列多行); 300+dpi; 紧裁白边
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.expanduser(os.environ.get(
    "B8GC_FIGDIR", os.path.join(os.environ.get("B8GC_WORK_ROOT", "./work"), "figures")))
os.makedirs(FIGDIR, exist_ok=True)
MM = 1 / 25.4
W1, W2 = 85 * MM, 170 * MM          # BMC 单栏/双栏 英寸
HMAX = 225 * MM                     # BMC 最大高

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 350,
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.titleweight": "bold", "figure.titlesize": 9,
})

STATE_COLORS = {"M0": "#6a5acd", "M1": "#2a9d8f", "M2": "#e76f51", "unrepro": "#cfcfcf"}
ROUTE_ORDER = ["primary", "ovarian_met", "ascites", "peritoneal_met"]
ROUTE_COLORS = {"primary": "#2a9d8f", "ovarian_met": "#8ab17d",
                "ascites": "#e9c46a", "peritoneal_met": "#e76f51"}
COHORT_CMAP = plt.cm.tab10


def panel(ax, letter):
    ax.text(-0.10, 1.08, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGDIR}/{name}.{ext}", bbox_inches="tight")
    print(f"saved {FIGDIR}/{name}.[pdf|png]")


def score_diff_axis(a, sigs):
    import scanpy as sc
    for s in ["M1", "M2"]:
        g = [x for x in sigs[s] if x in a.var_names]
        sc.tl.score_genes(a, g, score_name=f"s_{s}")
    return a.obs["s_M1"] - a.obs["s_M2"]


def morans_i(coords, x, k=6):
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    x = np.asarray(x, float); n = len(x)
    if n < k + 5:
        return np.nan
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords); idx = idx[:, 1:]
    xc = x - x.mean(); den = (xc ** 2).sum()
    if den == 0:
        return np.nan
    num = sum(xc[i] * xc[idx[i]].sum() for i in range(n))
    return (n / (n * k)) * (num / den)


def dotplot(ax, adata, genes, groupcol, groups, cmap="Reds"):
    """单细胞 dotplot: 点大小=表达比例, 颜色=z-scored 平均表达。"""
    import numpy as np
    import pandas as pd
    genes = [g for g in genes if g in adata.var_names]
    sub = adata[:, genes]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    df = pd.DataFrame(X, columns=genes); df["g"] = adata.obs[groupcol].astype(str).values
    mean = df.groupby("g")[genes].mean().reindex(groups)
    frac = df.groupby("g").apply(lambda d: (d[genes] > 0).mean()).reindex(groups)
    mz = (mean - mean.mean()) / (mean.std() + 1e-9)
    xs, ys, ss, cs = [], [], [], []
    for i, g in enumerate(groups):
        for j, gn in enumerate(genes):
            xs.append(j); ys.append(i); ss.append(frac.loc[g, gn] * 90 + 2); cs.append(mz.loc[g, gn])
    h = ax.scatter(xs, ys, s=ss, c=cs, cmap=cmap, vmin=-1, vmax=2, edgecolors="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=5)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups, fontsize=6)
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(-0.6, len(groups) - 0.4)
    return h


def feat_umap(ax, um, vals, title, cmap="viridis", vmax=None, s=1.2):
    """UMAP feature 散点(连续值)。vmax 可传分位数上限提对比。"""
    h = ax.scatter(um[:, 0], um[:, 1], c=vals, cmap=cmap, s=s, lw=0, vmax=vmax)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title)
    return h


def umap_labels(ax, um, labels, fs=6):
    """在每类质心直接标注名称(替代拥挤的 legend, 顶刊常用)。"""
    import numpy as np
    import pandas as pd
    lab = pd.Series(np.asarray(labels).astype(str))
    for c in lab.unique():
        mk = (lab == c).values
        cx, cy = np.median(um[mk, 0]), np.median(um[mk, 1])
        ax.text(cx, cy, c, fontsize=fs, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="0.6", lw=0.4, alpha=0.75))


def cbar(fig, h, ax, label=""):
    cb = fig.colorbar(h, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=5)
    if label:
        cb.set_label(label, fontsize=6)
    return cb
