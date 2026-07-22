#!/usr/bin/env python
"""s5_spatial.py — 空间验证（GSE251950 Visium）：M1/M2 是不是真空间组织域？

核心问题：分化态 M1/M2 在组织里是否形成空间分区（而非随机散布/解离假象）。
核心指标：分化轴 diff_axis(M1−M2) 的【空间自相关 Moran's I】——显著>0 = 分化在空间上成域。
辅助：M1 vs M2 在 spot 上的空间相关（负=互相排斥成不同域）。

⚠ 碰撞规避：原论文(Lee 2025 Gut)做 CCL2+成纤维×STAT3巨噬的免疫抑制 crosstalk；
  我们只问恶性细胞 M1/M2 空间分化，绝不碰它的 fibroblast/macrophage 故事。

签名来源：S3_5_states.h5ad（与 S3/S4 一致，去 MDK）。
输入  GSE251950 各 slide（自动搜 *filtered_feature_bc_matrix.h5 + tissue_positions）
Outputs under B8GC_WORK_ROOT/results/S5_spatial and B8GC_WORK_ROOT/qc.
"""
import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from sklearn.neighbors import NearestNeighbors

from sig_utils import EMT_INTRINSIC8

SPATIAL_ROOT = os.environ.get("B8GC_SPATIAL_ROOT", f"{C.HUB}/GC/GSE251950_raw")
OUTDIR = f"{C.RESULTS}/S5_spatial"
TABLE_S12 = f"{C.WORK}/figure_source_data/TableS12_spatial_morans_source.csv"
META = f"{C.RESULTS}/S3_5_states.h5ad"
STATES = ["M0", "M1", "M2"]
EMT_GENES = EMT_INTRINSIC8
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
MDK_AVOID = {"MDK"}


def derive_sigs():
    from sig_utils import frozen_sigs
    return frozen_sigs()


def load_slide(d):
    h5 = glob.glob(f"{d}/*filtered_feature_bc_matrix.h5")
    a = sc.read_10x_h5(h5[0]); a.var_names_make_unique()
    tp = (glob.glob(f"{d}/**/tissue_positions_list.csv", recursive=True)
          + glob.glob(f"{d}/**/tissue_positions.csv", recursive=True))
    raw = pd.read_csv(tp[0], header=None)
    if str(raw.iloc[0, 0]).lower() in ("barcode", "barcodes"):   # 有表头
        raw = pd.read_csv(tp[0])
    raw = raw.iloc[:, :6]
    raw.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    raw = raw.set_index("barcode")
    common = a.obs_names.intersection(raw.index)
    a = a[common].copy()
    a.obs["in_tissue"] = raw.loc[common, "in_tissue"].astype(int).values
    a.obsm["spatial"] = raw.loc[common, ["pxl_row", "pxl_col"]].astype(float).values
    a = a[a.obs["in_tissue"] == 1].copy()
    return a


def morans_i(coords, x, k=6):
    x = np.asarray(x, float)
    n = len(x)
    if n < k + 5:
        return np.nan
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    idx = idx[:, 1:]
    xc = x - x.mean()
    den = (xc ** 2).sum()
    if den == 0:
        return np.nan
    num = sum(xc[i] * xc[idx[i]].sum() for i in range(n))
    return (n / (n * k)) * (num / den)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(os.path.dirname(TABLE_S12), exist_ok=True)
    sigs = derive_sigs()
    slides = sorted({os.path.dirname(h) for h in
                     glob.glob(f"{SPATIAL_ROOT}/**/*filtered_feature_bc_matrix.h5", recursive=True)})
    print(f"找到 {len(slides)} 个 slide")

    rows = []
    for d in slides:
        name = os.path.basename(d)
        try:
            a = load_slide(d)
        except Exception as e:
            print(f"[skip] {name}: {e}"); continue
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        for s in STATES:
            gs = [g for g in sigs.get(s, []) if g in a.var_names]
            if gs:
                sc.tl.score_genes(a, gs, score_name=f"s_{s}")
        sc.tl.score_genes(a, [g for g in EMT_GENES if g in a.var_names], score_name="s_EMT")
        sc.tl.score_genes(a, [g for g in EPI if g in a.var_names], score_name="s_Epi")
        a.obs["diff_axis"] = a.obs.get("s_M1", 0) - a.obs.get("s_M2", 0)

        # 上皮高 spot（肿瘤区）
        epi_hi = a.obs["s_Epi"] > a.obs["s_Epi"].median()
        coords = a.obsm["spatial"]
        mi_all = morans_i(coords, a.obs["diff_axis"].values)
        mi_epi = morans_i(coords[epi_hi.values], a.obs["diff_axis"].values[epi_hi.values])
        mi_emt = morans_i(coords, a.obs["s_EMT"].values)
        # M1 vs M2 空间相关（负=排斥成不同域）
        m1m2_corr = np.corrcoef(a.obs.get("s_M1", pd.Series(0, index=a.obs.index)),
                                a.obs.get("s_M2", pd.Series(0, index=a.obs.index)))[0, 1]
        a.write(f"{OUTDIR}/{name}.h5ad")
        rows.append(dict(slide=name, n_spots=a.n_obs, n_epi_hi=int(epi_hi.sum()),
                         moransI_diffaxis=round(mi_all, 3), moransI_diffaxis_epi=round(mi_epi, 3),
                         moransI_EMT=round(mi_emt, 3), M1_M2_spatial_corr=round(m1m2_corr, 3)))
        print(f"  {name}: spots={a.n_obs} Moran'sI(diff)={mi_all:.3f} (epi {mi_epi:.3f}) "
              f"M1-M2corr={m1m2_corr:.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_S12, index=False)
    with open(f"{C.QC}/S5_summary.md", "w") as fh:
        fh.write("# S5 空间验证（Visium）：M1/M2 是否真空间组织域\n\n")
        fh.write(f"- {len(df)} 个 slide；Moran's I>0 且显著偏离 0 = 分化轴在空间成域\n\n")
        fh.write(df.to_markdown(index=False) + "\n\n")
        if len(df):
            fh.write(f"- **diff_axis Moran's I 中位数: {df['moransI_diffaxis'].median():.3f}**"
                     f"（上皮区 {df['moransI_diffaxis_epi'].median():.3f}）\n")
            fh.write(f"- EMT_intrinsic8 Moran's I 中位数: {df['moransI_EMT'].median():.3f}"
                     f"（intrinsic-8，与外部验证口径一致）\n")
            fh.write(f"- M1-M2 空间相关中位数: {df['M1_M2_spatial_corr'].median():.3f}"
                     f"（负=M1/M2 占不同空间域）\n")
    print("\n", df.to_string(index=False))
    print("\nS5 完成 → qc/S5_summary.md")


if __name__ == "__main__":
    main()
