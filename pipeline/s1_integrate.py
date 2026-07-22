#!/usr/bin/env python
"""s1_integrate.py — Stage 1：全细胞整合 + 大类注释（Harmony，CPU，256G 内存直接跑全量）。

Input: B8GC_WORK_ROOT/h5ad/<cohort>.h5ad (available cohorts; X is log-normalised).
Output: B8GC_WORK_ROOT/results/S1_integrated.h5ad
        obs: cohort/sample/route/leiden/cell_type ；obsm: X_pca / X_pca_harmony / X_umap ；raw=log归一化
      B8GC_WORK_ROOT/qc/S1_summary.md (cell_type by cohort summary)

Dependencies: harmonypy, leidenalg, and igraph (included in environment.yml).

设计：无 GPU → 用 Harmony(CPU) 做全细胞整合仅为“注释 + 联合可视化”；
      恶性细胞鉴定在 S2，复现性核心(MetaNeighbor 等)在 S3 的恶性子集上（不依赖整合）。
"""
import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import config as C

sc.settings.n_jobs = 16
sc.settings.verbosity = 1

OUT = f"{C.RESULTS}/S1_integrated.h5ad"
DROP_BLOOD = True          # route=blood(外周血)无恶性上皮 → 整合前剔除省算力；无该标签则自动忽略
N_HVG = 3000
N_PCS = 50
LEIDEN_RES = 1.0

# 大类 marker（注释用；Epithelial 含正常+恶性，S2 再分恶性）
MARKERS = {
    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "T_NK":        ["CD3D", "CD3E", "CD2", "NKG7", "GNLY"],
    "B":           ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma":      ["MZB1", "IGKC", "DERL3", "XBP1"],
    "Myeloid":     ["LYZ", "CD68", "CD14", "FCGR3A", "C1QA"],
    "Mast":        ["TPSAB1", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "ENG", "CLDN5"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "PDGFRB", "ACTA2"],
}


def load_all():
    ads = []
    for c in C.COHORTS:
        p = f"{C.H5AD}/{c}.h5ad"
        if not os.path.exists(p):
            print(f"  [skip] {c} (no h5ad)"); continue
        a = sc.read_h5ad(p)
        keep = [col for col in ("sample", "cohort", "route", "ln_status") if col in a.obs]
        a.obs = a.obs[keep].copy()
        a.layers.clear()                 # counts 各自 h5ad 留存；整合只用 log-norm X，省内存
        a.raw = None
        if "highly_variable" in a.var:   # 清掉 per-cohort HVG，后面统一重算
            a.var = a.var[[]].copy()
        ads.append(a); print(f"  loaded {c}: {a.n_obs} cells × {a.n_vars} genes")
    adata = sc.concat(ads, join="inner", index_unique="-")   # 基因取交集，干净对齐
    del ads; gc.collect()
    return adata


def main():
    os.makedirs(C.RESULTS, exist_ok=True)
    if os.path.exists(OUT):
        print(f"{OUT} 已存在 → 跳过 S1（要重跑先删它）"); return

    print("加载可用 cohort ...")
    adata = load_all()
    print(f"concat 后: {adata.shape}")

    if DROP_BLOOD and "route" in adata.obs and (adata.obs["route"] == "blood").any():
        n = int((adata.obs["route"] == "blood").sum())
        adata = adata[adata.obs["route"] != "blood"].copy()
        print(f"剔除 {n} 个外周血细胞 → {adata.n_obs}")

    adata.obs["cohort"] = adata.obs["cohort"].astype("category")

    # batch-aware HVG（log 数据上用 seurat flavor，避开 counts 整数要求）
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, batch_key="cohort", flavor="seurat")
    adata.raw = adata                                   # 存全基因 log-norm 供 score_genes

    # PCA（在 HVG 子集上 scale 后做，控制内存）
    hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=N_PCS)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    del hvg; gc.collect()

    # Harmony 整合（batch = cohort）
    print("Harmony 整合中（全细胞 CPU，约 0.5–1.5h）...")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, vars_use=["cohort"])
    z = np.asarray(ho.Z_corr)
    if z.shape != adata.obsm["X_pca"].shape:
        raise RuntimeError(f"Harmony output shape {z.shape} != PCA shape {adata.obsm['X_pca'].shape}")
    adata.obsm["X_pca_harmony"] = z

    # 邻居 + 聚类 + UMAP
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph",
                     n_iterations=2, directed=False)
    except Exception as e:
        print(f"  igraph leiden 失败({e})，回退 leidenalg"); sc.tl.leiden(adata, resolution=LEIDEN_RES)
    sc.tl.umap(adata)

    # 大类注释：每簇对各 marker set 打分，取最高分大类
    for ct, gs in MARKERS.items():
        gs = [g for g in gs if g in adata.raw.var_names]
        sc.tl.score_genes(adata, gs, score_name=f"sc_{ct}", use_raw=True)
    score_cols = [f"sc_{ct}" for ct in MARKERS]
    cl_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    cl2type = cl_scores.idxmax(axis=1).str.replace("sc_", "", regex=False)
    adata.obs["cell_type"] = adata.obs["leiden"].map(cl2type).astype("category")

    adata.write(OUT)

    ct_tab = pd.crosstab(adata.obs["cell_type"], adata.obs["cohort"])
    with open(f"{C.QC}/S1_summary.md", "w") as fh:
        fh.write("# S1 Integration Summary\n\n")
        fh.write(f"- cells: {adata.n_obs}; HVG-inner genes: {adata.n_vars}; "
                 f"leiden clusters: {adata.obs.leiden.nunique()}\n\n")
        fh.write("## cell_type × cohort（每个大类应跨多个队列出现=整合成功）\n\n")
        fh.write(ct_tab.to_markdown() + "\n\n")
        fh.write("## cell_type × route\n\n")
        fh.write(pd.crosstab(adata.obs["cell_type"], adata.obs["route"]).to_markdown() + "\n")
    print("S1 完成 →", OUT)
    print(ct_tab)


if __name__ == "__main__":
    main()
