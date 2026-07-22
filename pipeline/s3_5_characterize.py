#!/usr/bin/env python
"""s3_5_characterize.py — 刻画可复现的恶性 meta-cluster（回答"到底什么复现了"）。

为什么需要：S3v2 证明有 13 个簇跨队列复现，但 6 格预设程序太粗，把大多数塞进了
"Ribosomal" 垃圾桶。本步用真实 DE + GC 胃系谱系程序刻画这些可复现簇，得出真正的
"可复现恶性态分类"。

流程：
  1. 读 AUROC 矩阵(B8GC_MN_SCRATCH/auroc.csv) → 跨队列 AUROC>阈值建图 → 连通分量 = meta-cluster
  2. 保留跨 ≥MIN_COHORTS 队列的 meta-cluster = 可复现态
  3. 每个可复现 meta: 成分(队列/路径) + DE top markers + GC 胃系谱系打分 → 命名真实生物学
Outputs under B8GC_WORK_ROOT/results and B8GC_WORK_ROOT/qc.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

IN = f"{C.RESULTS}/S3v2_states.h5ad"
AUROC = os.path.join(os.environ.get("B8GC_MN_SCRATCH", os.path.join(C.WORK, "tmp", "mn")), "auroc.csv")
EDGE = 0.75
MIN_COHORTS = 3

# 扩展的 GC 恶性谱系程序（真实生物学，不止 6 格）
GC_LINEAGE = {
    "Pit_mucous_foveolar": ["MUC5AC", "TFF1", "GKN1", "GKN2", "FOXQ1"],
    "Neck_mucous": ["MUC6", "TFF2", "AQP5"],
    "Chief": ["PGA3", "PGC", "LIPF"],
    "Enteroendocrine": ["CHGA", "CHGB", "NEUROD1", "PCSK1N"],
    "Intestinal_metaplasia": ["CDX2", "MUC2", "REG4", "TFF3", "FABP1", "ITLN1"],
    "Proliferation": ["MKI67", "TOP2A", "CENPF", "UBE2C"],
    "EMT_pEMT": ["VIM", "ZEB1", "SPARC", "FN1", "TGFBI", "PDGFRA"],
    "Interferon": ["ISG15", "STAT1", "MX1", "IFIT3", "IFI6"],
    "Gastric_epithelial_general": ["EPCAM", "KRT8", "KRT18", "CLDN18", "KRT19"],
}


def build_metas(A, coh, nodes, edge):
    n = len(nodes)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if coh[i] != coh[j] and A[i, j] > edge:
                M[i, j] = M[j, i] = 1
    _, labels = connected_components(csr_matrix(M), directed=False)
    comp = pd.DataFrame({"node": nodes, "comp": labels, "cohort": coh})
    keep = [c for c, g in comp.groupby("comp") if g["cohort"].nunique() >= MIN_COHORTS]
    return comp, keep


def main():
    a = sc.read_h5ad(IN)
    a = a[a.obs["indep_cl"] != "NA"].copy()
    a.obs["node"] = a.obs["cohort"].astype(str) + "|" + a.obs["indep_cl"].astype(str)

    au = pd.read_csv(AUROC, index_col=0)
    nodes = [n for n in au.columns if n in set(a.obs["node"])]
    au = au.reindex(index=nodes, columns=nodes)
    A = np.nan_to_num((au.values + au.values.T) / 2.0)
    coh = np.array([n.split("|")[0] for n in nodes])

    # ---- 阈值扫描：更高 stringency 下大团会不会裂成多个跨队列态？----
    print("== AUROC 阈值扫描（≥3 队列的 meta 数 + 各 meta 的节点数）==")
    sweep_rows = []
    for edge in [0.75, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90]:
        comp, keep = build_metas(A, coh, nodes, edge)
        sizes = sorted([(comp["comp"] == c).sum() for c in keep], reverse=True)
        cohn = [comp[comp["comp"] == c]["cohort"].nunique() for c in keep]
        sweep_rows.append((edge, len(keep), sizes, sorted(cohn, reverse=True)))
        print(f"  edge={edge}: {len(keep)} 个≥3队列meta；节点数={sizes}；各覆盖队列数={sorted(cohn,reverse=True)}")

    # 选"最有结构"的阈值刻画：≥3队列meta 数最多者（并列取较高阈值）
    best = max(sweep_rows, key=lambda r: (r[1], r[0]))
    EDGE_USE = best[0]
    print(f"\n→ 用 edge={EDGE_USE} 刻画（该阈值下 {best[1]} 个跨队列 meta）")
    comp, keep = build_metas(A, coh, nodes, EDGE_USE)
    remap = {m: f"M{i}" for i, m in enumerate(sorted(set(keep)))}
    node2meta = dict(zip(comp["node"], comp["comp"]))
    a.obs["meta"] = (a.obs["node"].map(node2meta)
                     .map(lambda x: remap.get(x, "unrepro")).astype("category"))
    n_meta = len(keep)

    # 谱系打分
    for p, gs in GC_LINEAGE.items():
        gs = [g for g in gs if g in a.var_names]
        if gs:
            sc.tl.score_genes(a, gs, score_name=f"L_{p}", use_raw=False)
    lcols = [f"L_{p}" for p in GC_LINEAGE if f"L_{p}" in a.obs]

    # DE：可复现 meta 之间
    sub = a[a.obs["meta"] != "unrepro"].copy()
    de_markers = {}
    if sub.obs["meta"].nunique() >= 2:
        sc.tl.rank_genes_groups(sub, "meta", method="wilcoxon", n_genes=25)
        for m in sub.obs["meta"].cat.categories:
            if m == "unrepro":
                continue
            de_markers[m] = list(sub.uns["rank_genes_groups"]["names"][m][:15])

    a.write(f"{C.RESULTS}/S3_5_states.h5ad")
    with open(f"{C.QC}/S3_5_summary.md", "w") as fh:
        fh.write("# S3.5 可复现恶性 meta-cluster 刻画\n\n")
        fh.write(f"- 可复现 meta-cluster(跨≥{MIN_COHORTS}队列): **{n_meta} 个**\n")
        fh.write(f"- AUROC 建图阈值 {EDGE}；总细胞 {a.n_obs}\n\n")
        for m in remap.values():
            cells = a[a.obs["meta"] == m]
            fh.write(f"## {m}（{cells.n_obs} cells）\n\n")
            fh.write(f"- 队列: {dict(cells.obs.cohort.value_counts()[cells.obs.cohort.value_counts()>0])}\n")
            fh.write(f"- 路径: {dict(cells.obs.route.value_counts()[cells.obs.route.value_counts()>0])}\n")
            ls = cells.obs[lcols].mean().sort_values(ascending=False)
            top_prog = ls.index[0].replace("L_", "")
            fh.write(f"- 最高谱系程序: **{top_prog}** (score {ls.iloc[0]:.3f})；"
                     f"次高 {ls.index[1].replace('L_','')} ({ls.iloc[1]:.3f})\n")
            fh.write(f"- DE top markers: {', '.join(de_markers.get(m, [])[:12])}\n\n")
        fh.write("## 谱系程序 × meta（每行一个可复现态的平均谱系打分）\n\n")
        prof = a[a.obs.meta != "unrepro"].obs.groupby("meta", observed=True)[lcols].mean()
        prof.columns = [c.replace("L_", "") for c in prof.columns]
        fh.write(prof.round(3).to_markdown() + "\n")
    print("S3.5 完成。谱系 profile:")
    print(a[a.obs.meta != "unrepro"].obs.groupby("meta", observed=True)[lcols].mean().round(3))


if __name__ == "__main__":
    main()
