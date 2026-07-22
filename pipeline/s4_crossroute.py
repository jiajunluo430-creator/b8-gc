#!/usr/bin/env python
"""s4_crossroute.py — 分化轴的跨路径持续 + EMT 不随转移富集（论文核心正面图 + punchline）。

骨架（来自 S3.6 稳健性结论）：
  可复现分化轴 = 黏液分化态 M1（TFF1/TFF2/MUC5AC）↔ 低分化/分泌态 M2（PCSK1N/WFDC2/CLU）；
  M0（小凹-代谢）为次要 facet。问：这条轴在原发→卵巢→腹膜转移里是否保留？EMT 是否随转移富集？

关键诚实处理：
  - route 与 cohort 混杂 → 跨路径组成【主分析在 HRA004702 内部】（唯一多路径恶性队列），
    其余队列只作为"态本身可复现"的锚（S3 已证）。
  - M2 签名含 MDK（HRA 的 ER-MDK-LRP1 机制基因）→ 本步不以 MDK 描述 M2，规避碰撞。
  - EMT 用预设基因集打分，证其不构成离散可复现态、且不随转移系统性升高。

输入  S3_5_states.h5ad（meta 标签，导签名）, S2_malignant.h5ad（全部恶性）
Outputs under B8GC_WORK_ROOT/results and B8GC_WORK_ROOT/qc.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import config as C

IN_META = f"{C.RESULTS}/S3_5_states.h5ad"
MAL = f"{C.RESULTS}/S2_malignant.h5ad"
CARC = ["primary", "ovarian_met", "peritoneal_met", "ascites", "peritoneal_lavage"]
STATES = ["M0", "M1", "M2"]
EMT_GENES = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2",
             "FN1", "SPARC", "TGFBI", "ITGA5", "TIMP1", "MMP2", "LOXL2"]
MDK_AVOID = {"MDK"}    # 不进 M2 描述（碰撞规避）


def derive_sigs(meta):
    from sig_utils import frozen_sigs
    sigs = frozen_sigs()
    have = [m for m in STATES if m in sigs]
    return {m: [x for x in sigs[m] if x not in MDK_AVOID] for m in have}, have


def main():
    meta = sc.read_h5ad(IN_META)
    sigs, have = derive_sigs(meta)

    mal = sc.read_h5ad(MAL)
    mal = mal[mal.obs["route"].isin(CARC)].copy()
    scols = []
    for m in have:
        gs = [g for g in sigs[m] if g in mal.var_names]
        sc.tl.score_genes(mal, gs, score_name=f"s_{m}"); scols.append(f"s_{m}")
    egs = [g for g in EMT_GENES if g in mal.var_names]
    sc.tl.score_genes(mal, egs, score_name="s_EMT")
    mal.obs["state"] = mal.obs[scols].idxmax(axis=1).str.replace("s_", "", regex=False)

    # 分化轴标量：M1(分化) - M2(低分化)
    mal.obs["diff_axis"] = mal.obs["s_M1"] - mal.obs["s_M2"]

    L = ["# S4 分化轴跨路径持续 + EMT 不随转移富集\n"]

    # ---- 主分析：HRA004702 内部（唯一多路径恶性队列）----
    hra = mal[mal.obs["cohort"] == "HRA004702"]
    L += ["## 主分析 — HRA004702 内部（控住 cohort，看纯 route 效应）\n",
          "### 态组成 × route（argmax 占比）\n",
          pd.crosstab(hra.obs["route"], hra.obs["state"], normalize="index").round(3).to_markdown(), "\n",
          "### 平均签名分 × route（连续，更稳）\n",
          hra.obs.groupby("route", observed=True)[scols + ["s_EMT", "diff_axis"]].mean().round(3).to_markdown(), "\n"]

    # ---- 跨队列汇总（注明 cohort 混杂）----
    L += ["## 跨队列汇总 route 组成（⚠cohort 混杂，仅辅助）\n",
          pd.crosstab(mal.obs["route"], mal.obs["state"], normalize="index").round(3).to_markdown(), "\n"]

    # ---- punchline: EMT 不构成转移驱动轴、且不随转移富集 ----
    emt_by_state = mal.obs.groupby("state", observed=True)["s_EMT"].mean().round(3)
    emt_by_route_hra = hra.obs.groupby("route", observed=True)["s_EMT"].agg(["mean", "median"]).round(3)
    # EMT 高细胞(top 10%)在各 route 的占比 —— 若不随转移升高，则 EMT 非转移可复现态
    thr = mal.obs["s_EMT"].quantile(0.9)
    mal.obs["EMT_hi"] = mal.obs["s_EMT"] > thr
    emt_hi_route = mal.obs.groupby("route", observed=True)["EMT_hi"].mean().round(3)
    L += ["## ★Punchline — EMT 不随转移富集（非转移驱动轴）\n",
          f"- EMT 分 × 态（无任一态以 EMT 为主 → EMT 不构成离散态）:\n\n{emt_by_state.to_markdown()}\n",
          f"- HRA 内 EMT 分 × route（不随转移单调升高 → EMT 非转移驱动轴）:\n\n{emt_by_route_hra.to_markdown()}\n",
          f"- EMT-high(top10%) 细胞各 route 占比:\n\n{emt_hi_route.to_markdown()}\n"]

    # ---- 分化轴跨路径持续小结 ----
    L += ["## 分化轴跨路径持续（M1−M2 = diff_axis；正=偏分化）\n",
          hra.obs.groupby("route", observed=True)["diff_axis"].agg(["mean", "median"]).round(3).to_markdown(), "\n",
          "解读提示：若 M1/M2 在转移 route 仍都占可观比例 → 分化轴持续到转移；"
          "若 diff_axis 在转移端不崩 → 转移未全面去分化。\n"]

    mal.write(f"{C.RESULTS}/S4_scored.h5ad")
    with open(f"{C.QC}/S4_summary.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nS4 完成 → qc/S4_summary.md, results/S4_scored.h5ad")


if __name__ == "__main__":
    main()
