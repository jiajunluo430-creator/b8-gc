#!/usr/bin/env python
"""s3_6_robustness.py — 可复现态稳健性自检（进 S4 前的把关）。

回答两个问题：
  Q1. M0 vs M1 是不是两个真态，还是该合并成一个黏液/小凹态？
      → 直接 DE，数强差异基因；<阈值则建议合并。
  Q2. 这 3 个态的签名，在【每个队列各自】是否都能稳定标出一群细胞？
      → 从 meta 细胞导出签名 → 在各队列【全部】恶性细胞上独立打分 → 每队列是否都有该态。
      这才是"跨队列可复现"的硬证据（不依赖之前的聚类/阈值）。

输入  S3_5_states.h5ad（带 meta 标签）, S2_malignant.h5ad（全部恶性，用于跨队列打分）
Output: B8GC_WORK_ROOT/qc/S3_6_robustness.md
"""
import numpy as np
import pandas as pd
import scanpy as sc
import config as C

IN = f"{C.RESULTS}/S3_5_states.h5ad"
MAL = f"{C.RESULTS}/S2_malignant.h5ad"
CARC = ["primary", "ovarian_met", "peritoneal_met", "ascites", "peritoneal_lavage"]
METAS = ["M0", "M1", "M2"]
DE_FDR, DE_LFC = 0.05, 1.0
MIN_FRAC = 0.02          # 某态在某队列占比 ≥2% 才算"该队列存在此态"


def main():
    a = sc.read_h5ad(IN)
    have = [m for m in METAS if (a.obs["meta"] == m).any()]
    lines = ["# S3.6 可复现态稳健性自检\n"]

    # ---- Q1: M0 vs M1 是否该合并 ----
    if "M0" in have and "M1" in have:
        m01 = a[a.obs["meta"].isin(["M0", "M1"])].copy()
        m01.obs["meta"] = m01.obs["meta"].astype(str)
        sc.tl.rank_genes_groups(m01, "meta", method="wilcoxon")
        de = sc.get.rank_genes_groups_df(m01, None)
        strong = de[(de["pvals_adj"] < DE_FDR) & (de["logfoldchanges"].abs() > DE_LFC)]
        n0 = int((strong["group"] == "M0").sum())
        n1 = int((strong["group"] == "M1").sum())
        top0 = list(strong[strong.group == "M0"].sort_values("logfoldchanges", ascending=False)["names"][:10])
        top1 = list(strong[strong.group == "M1"].sort_values("logfoldchanges", ascending=False)["names"][:10])
        merge = (n0 < 10) or (n1 < 10)
        lines += [
            "## Q1: M0 vs M1 区分度\n",
            f"- 强差异基因(FDR<{DE_FDR}, |log2FC|>{DE_LFC}): M0 侧 {n0} 个, M1 侧 {n1} 个",
            f"- M0 top: {', '.join(top0)}",
            f"- M1 top: {', '.join(top1)}",
            f"- **结论: {'两侧都太弱 → 建议合并为 1 个黏液/小凹态' if merge else '两侧都有足够特异基因 → 保留为 2 个态'}**\n",
        ]
    else:
        lines += ["## Q1: M0/M1 至少一个不存在，跳过合并判断\n"]

    # ---- 导出各态签名（meta vs rest）----
    a.obs["meta"] = a.obs["meta"].astype(str)
    sc.tl.rank_genes_groups(a, "meta", groups=have, reference="rest", method="wilcoxon")
    sigs = {}
    for m in have:
        d = sc.get.rank_genes_groups_df(a, m)
        sigs[m] = list(d[(d["pvals_adj"] < DE_FDR) & (d["logfoldchanges"] > DE_LFC)]["names"][:30])

    # ---- Q2: 签名在各队列独立打分 → 跨队列存在性 ----
    mal = sc.read_h5ad(MAL)
    mal = mal[mal.obs["route"].isin(CARC)].copy()
    scol = []
    for m in have:
        gs = [g for g in sigs[m] if g in mal.var_names]
        if gs:
            sc.tl.score_genes(mal, gs, score_name=f"sig_{m}"); scol.append(f"sig_{m}")
    mal.obs["assigned"] = mal.obs[scol].idxmax(axis=1).str.replace("sig_", "", regex=False)
    frac = pd.crosstab(mal.obs["cohort"], mal.obs["assigned"], normalize="index").round(3)
    corr = mal.obs[scol].corr().round(2)
    # 每个态在多少队列里占比≥MIN_FRAC
    present = {m: int((frac.get(m, pd.Series(dtype=float)) >= MIN_FRAC).sum()) for m in have}

    lines += [
        "## Q2: 签名跨队列存在性（每队列恶性细胞按签名 argmax 归属的占比）\n",
        frac.to_markdown(), "\n",
        f"- 各态出现在 ≥{int(MIN_FRAC*100)}% 占比的队列数: {present}",
        f"- 签名两两相关(低=互相独立、是真不同态):\n\n{corr.to_markdown()}\n",
    ]
    repro = [m for m in have if present.get(m, 0) >= 3]
    lines += [
        f"## 总判决\n",
        f"- 跨 ≥3 队列稳定存在的态: **{repro}**（共 {len(repro)} 个）",
        f"- 每态签名(top markers):",
    ]
    for m in have:
        lines.append(f"  - {m}: {', '.join(sigs[m][:15])}")

    with open(f"{C.QC}/S3_6_robustness.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\nS3.6 完成 → qc/S3_6_robustness.md")


if __name__ == "__main__":
    main()
