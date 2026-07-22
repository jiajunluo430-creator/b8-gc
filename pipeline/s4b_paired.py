#!/usr/bin/env python
"""s4b_paired.py — 同病人配对：分化轴跨路径变化是"选择"还是"转化"？

S4 发现卵巢转移保留分化、腹膜转移去分化，但那是跨样本比，可能只是病人选择偏差。
HRA004702 是配对设计（同病人 P01_PT / P01_OM / P01_PM / P01_AS），用【同病人内】
原发 vs 转移的 diff_axis 配对检验，才能证明是转移导致的转化，而非选择。

输入  S4_scored.h5ad（含 diff_axis / s_M1 / s_M2 / state / sample / route / cohort）
Output: B8GC_WORK_ROOT/qc/S4b_paired.md
"""
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from scipy.stats import wilcoxon

IN = f"{C.RESULTS}/S4_scored.h5ad"
METS = ["ovarian_met", "peritoneal_met", "ascites"]


def main():
    a = sc.read_h5ad(IN)
    hra = a[a.obs["cohort"] == "HRA004702"].copy()
    # 病人 = sample 去掉 route 后缀（P01_OM → P01）
    hra.obs["patient"] = hra.obs["sample"].astype(str).str.rsplit("_", n=1).str[0]
    print("HRA samples:", sorted(hra.obs["sample"].unique())[:12], "...")
    print("HRA patients:", sorted(hra.obs["patient"].unique()))

    # 每病人 × route 的平均 diff_axis 与 M2 占比
    pr_diff = hra.obs.groupby(["patient", "route"], observed=True)["diff_axis"].mean().unstack()
    hra.obs["is_M2"] = (hra.obs["state"] == "M2").astype(int)
    pr_m2 = hra.obs.groupby(["patient", "route"], observed=True)["is_M2"].mean().unstack()

    L = ["# S4b 同病人配对：分化轴 选择 vs 转化\n"]
    L += [f"- 病人数: {hra.obs['patient'].nunique()}；samples: {hra.obs['sample'].nunique()}\n"]
    L += ["## 每病人 × route 平均 diff_axis（M1−M2，正=分化）\n",
          pr_diff.round(3).to_markdown(), "\n"]

    L += ["## 同病人配对检验（原发 vs 各转移，diff_axis）\n"]
    for met in METS:
        if "primary" not in pr_diff.columns or met not in pr_diff.columns:
            L.append(f"- primary vs {met}: 缺列，跳过"); continue
        sub = pr_diff[["primary", met]].dropna()
        n = len(sub)
        if n < 3:
            L.append(f"- primary vs {met}: 仅 {n} 个配对病人，样本太少（仅描述）"
                     f" 均值 primary={sub['primary'].mean():.3f} {met}={sub[met].mean():.3f}")
            continue
        try:
            stat, p = wilcoxon(sub["primary"], sub[met])
        except Exception:
            p = np.nan
        delta = (sub[met] - sub["primary"]).mean()
        direction = "去分化" if delta < 0 else "更分化"
        L.append(f"- **primary vs {met}**: n={n} 配对病人；"
                 f"原发均值 {sub['primary'].mean():.3f} → {met} 均值 {sub[met].mean():.3f}；"
                 f"Δ={delta:+.3f}（{direction}）；Wilcoxon p={p:.4f}")

    L += ["\n## 每病人 × route 的 M2(低分化)占比\n", pr_m2.round(3).to_markdown(), "\n"]
    L += ["\n解读：若同病人 primary→peritoneal_met 的 diff_axis 显著下降(p<0.05, Δ<0) → "
          "去分化是【转化】而非选择，punchline 成立；若 primary→ovarian_met 无显著下降 → "
          "卵巢转移确实【保留】分化。配对病人少时 p 仅供参考，以 Δ 方向 + 个体一致性为主。"]

    with open(f"{C.QC}/S4b_paired.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nS4b 完成 → qc/S4b_paired.md")


if __name__ == "__main__":
    main()
