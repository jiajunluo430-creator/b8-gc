#!/usr/bin/env python
"""s3_reproducibility_v2.py — Stage 3 修正版（去循环论证 + 剔技术伪迹）。★G2 生死闸★

为什么重写 v1：v1 把态定义在 Harmony 联合空间，再问"跨队列一致吗"——整合已强制对齐，
AUROC 必然≈1（循环）。且一半态是 Stress/OxPhos 技术伪迹。v2 三处修正：

  1. 去循环：每个队列【独立】聚类 → 独立 state 标签；MetaNeighbor 跨队列比对独立簇，
     高 AUROC = 同一个态在不同队列里被【各自独立】发现 = 真复现。整合不参与判据。
  2. 只用癌组织：primary + 各转移路径；precancer / adjacent_normal 排除出核心
     （precancer 另存做进展参照；adjacent_normal 的"恶性"多为 CopyKAT 假阳）。
  3. 标技术 vs 生物：Cycling/Stress/OxPhos/Ribosomal=技术或trivial；只有
     EMT_pEMT/Glandular_GI/Interferon 等生物程序的可复现态才算数。

G2 v2 判决：≥1 个【生物学】meta-program 在 ≥3 队列里被独立发现且互配 AUROC>0.85 → PASS。

依赖: pip install pymn scanorama
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from sklearn.metrics import adjusted_rand_score

sc.settings.n_jobs = 16
IN = f"{C.RESULTS}/S2_malignant.h5ad"
OUT = f"{C.RESULTS}/S3v2_states.h5ad"
CARCINOMA = ["primary", "ovarian_met", "peritoneal_met", "ascites", "peritoneal_lavage"]
PERCOHORT_RES = 0.5
AUROC_MATCH = 0.85          # 跨队列独立簇互配阈值
MIN_COHORTS = 3
MIN_CL = 50                 # 队列内独立簇最小细胞数

TECH_PROGRAMS = {"Cycling", "Stress", "OxPhos", "Ribosomal"}   # 技术/trivial
PROG_MARKERS = {
    "Cycling": ["MKI67", "TOP2A", "CENPF", "UBE2C"],
    "Stress": ["HSPA1A", "HSPA1B", "JUN", "FOS", "DNAJB1", "HSPB1"],
    "OxPhos": ["NDUFA4", "COX5B", "ATP5F1E", "COX6C"],
    "Ribosomal": ["RPS6", "RPL13", "RPS18", "RPL10"],
    "EMT_pEMT": ["VIM", "ZEB1", "SPARC", "TIMP1", "FN1", "TGFBI"],
    "Glandular_GI": ["TFF1", "TFF3", "MUC5AC", "REG4", "LGR5", "MUC2"],
    "Interferon": ["ISG15", "IFI6", "STAT1", "MX1", "IFIT3"],
}


def is_tech_gene(g):
    """核糖体/线粒体/血红蛋白/热休克/即刻早期/核 lncRNA → 技术轴，聚类与 MetaNeighbor 都排除。"""
    g = str(g).upper()
    if g.startswith(("RPS", "RPL", "MRPS", "MRPL", "MT-", "HBA", "HBB",
                     "HSPA", "HSPB", "HSPH", "HSPD", "HSPE", "DNAJ")):
        return True
    return g in {"FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "DUSP1",
                 "HSP90AA1", "HSP90AB1", "MALAT1", "NEAT1", "BAG3", "JCHAIN", "XIST"}


def per_cohort_independent_clusters(a):
    """每个队列各自聚类（绝不联合整合）→ obs['indep_cl']（队列内编号）。"""
    a.obs["indep_cl"] = "NA"
    for coh in a.obs["cohort"].cat.categories:
        idx = a.obs["cohort"] == coh
        if idx.sum() < MIN_CL * 2:
            continue
        sub = a[idx].copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
        # 关键：从聚类特征里剔除技术轴基因（核糖体/线粒体/应激…），否则技术梯度主导聚类
        tech = sub.var_names.to_series().map(is_tech_gene).values
        sub.var.loc[tech, "highly_variable"] = False
        h = sub[:, sub.var.highly_variable].copy()
        sc.pp.scale(h, max_value=10)
        sc.tl.pca(h, n_comps=30)
        sc.pp.neighbors(h, use_rep="X_pca", n_neighbors=15)
        try:
            sc.tl.leiden(h, resolution=PERCOHORT_RES, flavor="igraph",
                         n_iterations=2, directed=False)
        except Exception:
            sc.tl.leiden(h, resolution=PERCOHORT_RES)
        # 丢弃过小的独立簇
        vc = h.obs["leiden"].value_counts()
        keep = set(vc[vc >= MIN_CL].index)
        lab = h.obs["leiden"].astype(str).where(h.obs["leiden"].isin(keep), other="NA")
        a.obs.loc[idx, "indep_cl"] = lab.values
    return a


def annotate_programs(a):
    for p, gs in PROG_MARKERS.items():
        gs = [g for g in gs if g in a.var_names]
        if gs:
            sc.tl.score_genes(a, gs, score_name=f"p_{p}", use_raw=False)
    pc = [f"p_{p}" for p in PROG_MARKERS if f"p_{p}" in a.obs]
    a.obs["node"] = a.obs["cohort"].astype(str) + "|" + a.obs["indep_cl"].astype(str)
    prog = a.obs.groupby("node", observed=True)[pc].mean().idxmax(axis=1).str.replace("p_", "")
    return prog            # node -> program


def metaneighbor_independent(a):
    """导出独立簇表达+标签 → 调 canonical R MetaNeighbor → 读回 AUROC 矩阵。"""
    import subprocess
    import scipy.io as sio
    import scipy.sparse as sp
    sub = a[a.obs["indep_cl"] != "NA"].copy()
    # 同样从 MetaNeighbor 的基因集里剔除技术轴，让 variableGenes 只挑生物学基因
    keep = ~sub.var_names.to_series().map(is_tech_gene).values
    sub = sub[:, keep].copy()
    print(f"[MetaNeighbor] 剔除技术基因后剩 {sub.n_vars} 基因")
    mn_dir = os.environ.get("B8GC_MN_SCRATCH", os.path.join(C.WORK, "tmp", "mn"))
    os.makedirs(mn_dir, exist_ok=True)
    X = sub.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    sio.mmwrite(os.path.join(mn_dir, "expr.mtx"), X.T.tocoo())    # genes × cells
    pd.Series(sub.var_names).to_csv(os.path.join(mn_dir, "genes.txt"), index=False, header=False)
    sub.obs[["cohort", "indep_cl"]].astype(str).to_csv(os.path.join(mn_dir, "meta.csv"), index=False)
    rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_metaneighbor.R")
    rscript_bin = os.environ.get("B8GC_RSCRIPT", "Rscript")
    print(f"[MetaNeighbor] 调 R: {rscript}")
    r = subprocess.run([rscript_bin, rscript, mn_dir], capture_output=True, text=True)
    print(r.stdout)
    auroc_path = os.path.join(mn_dir, "auroc.csv")
    if r.returncode != 0 or not os.path.exists(auroc_path):
        raise RuntimeError(f"R MetaNeighbor 失败：\n{r.stderr}\n"
                           f"（缺包则先 R -e 'BiocManager::install(\"MetaNeighbor\")'）")
    return pd.read_csv(auroc_path, index_col=0)


def reciprocal_hits(auroc):
    """对每个节点找跨队列最佳匹配；统计每个节点在多少其它队列里有 AUROC>阈值 的匹配。"""
    nodes = list(auroc.columns)
    coh = {n: n.split("|")[0] for n in nodes}
    rows = []
    for n in nodes:
        others = {}
        for m in nodes:
            if coh[m] == coh[n]:
                continue
            others.setdefault(coh[m], []).append(auroc.loc[n, m])
        best_per_coh = {c: max(v) for c, v in others.items()}
        hit_cohorts = [c for c, v in best_per_coh.items() if v > AUROC_MATCH]
        rows.append(dict(node=n, cohort=coh[n],
                         n_hit_cohorts=len(hit_cohorts),
                         best_match_auroc=max(best_per_coh.values()) if best_per_coh else np.nan,
                         median_cross_auroc=float(np.median(list(best_per_coh.values())))
                         if best_per_coh else np.nan))
    return pd.DataFrame(rows)


def main():
    a = sc.read_h5ad(IN)
    a.obs["cohort"] = a.obs["cohort"].astype("category")
    print("S2 恶性总数:", a.n_obs, dict(a.obs["route"].value_counts()))

    # precancer 单独存（进展参照），核心只留癌组织
    a[a.obs["route"] == "precancer"].copy().write(f"{C.RESULTS}/S3v2_precancer.h5ad")
    car = a[a.obs["route"].isin(CARCINOMA)].copy()
    print(f"核心癌组织恶性: {car.n_obs}（排除 precancer/adjacent_normal）",
          dict(car.obs["route"].value_counts()))

    car = per_cohort_independent_clusters(car)
    prog = annotate_programs(car)
    auroc = metaneighbor_independent(car)
    hits = reciprocal_hits(auroc)
    hits["program"] = hits["node"].map(prog)
    hits["is_technical"] = hits["program"].isin(TECH_PROGRAMS)

    # 可复现 = 在≥(MIN_COHORTS-1)个其它队列里有匹配
    hits["reproducible"] = hits["n_hit_cohorts"] >= (MIN_COHORTS - 1)
    bio_repro = hits[(hits.reproducible) & (~hits.is_technical)]
    bio_programs = sorted(bio_repro["program"].unique())
    verdict = "PASS" if len(bio_programs) >= 1 else "FAIL/审查"

    # 双整合 ARI（在核心癌组织上，粗分辨）——稳健性佐证
    ari = np.nan
    try:
        import scanorama, anndata as ad
        sp_ = [car[car.obs["cohort"] == c].copy() for c in car.obs["cohort"].cat.categories
               if (car.obs["cohort"] == c).sum() > 0]
        scanorama.integrate_scanpy(sp_)
        mg = ad.concat(sp_, join="inner")[car.obs_names].copy()
        sc.pp.neighbors(mg, use_rep="X_scanorama")
        sc.tl.leiden(mg, resolution=0.3, key_added="scan")
        # Harmony 粗聚类
        sc.pp.highly_variable_genes(car, n_top_genes=2000, batch_key="cohort", flavor="seurat")
        hh = car[:, car.var.highly_variable].copy(); sc.pp.scale(hh, max_value=10)
        sc.tl.pca(hh, n_comps=30); car.obsm["X_pca"] = hh.obsm["X_pca"]
        sc.external.pp.harmony_integrate(car, "cohort", basis="X_pca", adjusted_basis="X_h")
        sc.pp.neighbors(car, use_rep="X_h"); sc.tl.leiden(car, resolution=0.3, key_added="harm")
        ari = adjusted_rand_score(car.obs["harm"], mg.obs["scan"])
    except Exception as e:
        print(f"[ARI] 跳过: {e}")

    car.write(OUT)
    auroc.to_csv(f"{C.RESULTS}/auroc.csv")
    auroc.to_csv(f"{C.RESULTS}/S3v2_metaneighbor_auroc.csv")
    hits.to_csv(f"{C.RESULTS}/S3v2_reproducibility_hits.csv", index=False)
    with open(f"{C.QC}/S3v2_summary.md", "w") as fh:
        fh.write("# S3 v2 Reproducibility (non-circular) — Gate G2\n\n")
        fh.write(f"核心癌组织恶性 {car.n_obs}，{car.obs['cohort'].nunique()} 队列；"
                 f"独立簇节点 {auroc.shape[0]} 个\n\n")
        fh.write(f"**双整合 ARI(粗分辨) = {ari}**（>0.6 为稳健）\n\n")
        fh.write("## 独立簇跨队列复现（每个节点=某队列独立聚出的簇）\n\n")
        fh.write(hits.sort_values(["reproducible", "n_hit_cohorts", "median_cross_auroc"],
                                  ascending=False).to_markdown(index=False) + "\n\n")
        fh.write(f"## G2 v2 判决: **{verdict}**\n\n")
        fh.write(f"- 可复现【生物学】program（≥{MIN_COHORTS}队列独立发现, 非技术）: "
                 f"{bio_programs if bio_programs else '无'}\n")
        fh.write(f"- 被判技术伪迹而剔除的可复现节点: "
                 f"{sorted(hits[(hits.reproducible)&(hits.is_technical)]['program'].unique())}\n")
    print(f"\n=== G2 v2 判决: {verdict} ===  生物学可复现 program: {bio_programs}")
    print(hits.sort_values(["reproducible", "n_hit_cohorts"], ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
