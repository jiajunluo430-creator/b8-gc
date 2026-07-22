#!/usr/bin/env python
"""s6b_grade_adjust.py — 预后这块的最后两个测试（不 fishing，出什么写什么）。

① diff_axis vs 病理分化 Grade（TCGA）：计算分化轴是否吻合病理等级 → 签名的正交验证。
② stage 校正的多变量 Cox：校正分期/年龄后 diff_axis 还有没有独立预后。
两者出什么报什么，不再调 cutoff、不再逐个试 M0/M1/M2。

Output: B8GC_WORK_ROOT/qc/S6b_grade_adjust.md
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from scipy.stats import spearmanr, kruskal
from lifelines import CoxPHFitter

META = f"{C.RESULTS}/S3_5_states.h5ad"
EXPR = f"{C.HUB}/GC/TCGA_STAD_Exp.txt"
CLIN = f"{C.HUB}/GC/TCGA_STAD_clinical.txt"
TIME = f"{C.HUB}/GC/TCGA_STAD_time.txt"
MDK_AVOID = {"MDK"}


def derive_sigs():
    from sig_utils import frozen_sigs
    return frozen_sigs()


def grade_ord(g):
    return {"G1": 1, "G2": 2, "G3": 3, "G4": 4}.get(str(g).strip().upper(), np.nan)


def stage_ord(s):
    s = str(s).upper().replace("STAGE", "").strip()
    for k, v in [("IV", 4), ("III", 3), ("II", 2), ("I", 1)]:
        if s.startswith(k):
            return v
    return np.nan


def age_hi(a):
    s = str(a).strip().upper()
    if s == ">65":
        return 1
    if s == "<=65":
        return 0
    try:
        return int(float(s) > 65)
    except Exception:
        return np.nan


def main():
    sigs = derive_sigs()
    e = pd.read_csv(EXPR, sep=",").drop(columns=["gene_id", "gene_type"], errors="ignore").set_index("gene_name")
    e = e.apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    e = e.loc[:, [c for c in e.columns if str(c).endswith("-01")]]
    e.columns = ["-".join(str(c).split("-")[:3]) for c in e.columns]
    e = e.loc[:, ~pd.Index(e.columns).duplicated()]
    if np.nanmax(e.values) > 100:
        e = np.log1p(e)

    def score(genes):
        g = [x for x in genes if x in e.index]
        z = e.loc[g]; z = z.sub(z.mean(1), axis=0).div(z.std(1) + 1e-9, axis=0)
        return z.mean(0)
    diff = (score(sigs["M1"]) - score(sigs["M2"])).rename("diff_axis")
    df = diff.reset_index(); df.columns = ["sample", "diff_axis"]

    clin = pd.read_csv(CLIN, sep="\t").rename(columns={"Id": "sample"})
    clin["sample"] = clin["sample"].astype(str)
    clin["grade_ord"] = clin["Grade"].map(grade_ord)
    clin["stage_ord"] = clin["Stage"].map(stage_ord)
    clin["age_hi"] = clin["Age"].map(age_hi)
    tm = pd.read_csv(TIME, sep="\t").rename(columns={"Id": "sample", "futime": "time", "fustat": "event"})
    tm["sample"] = tm["sample"].astype(str)

    d = df.merge(clin[["sample", "grade_ord", "stage_ord", "age_hi"]], on="sample").merge(
        tm[["sample", "time", "event"]], on="sample")
    for c in ["time", "event"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    L = ["# S6b 分化签名的 Grade 验证 + stage 校正 Cox（TCGA-STAD）\n"]

    # ① diff_axis vs Grade
    g = d.dropna(subset=["grade_ord", "diff_axis"])
    rho, pg = spearmanr(g["grade_ord"], g["diff_axis"])
    by = g.groupby("grade_ord")["diff_axis"].agg(["mean", "count"]).round(3)
    try:
        kw = kruskal(*[grp["diff_axis"].values for _, grp in g.groupby("grade_ord")])
        kwp = kw.pvalue
    except Exception:
        kwp = np.nan
    L += ["## ① diff_axis vs 病理分化 Grade（验证签名生物学意义）\n",
          f"- Spearman rho={rho:.3f}, p={pg:.4f}（**负=高级别低分化→低 diff_axis，吻合预期**）",
          f"- Kruskal-Wallis p={kwp:.4f}",
          f"- 各 Grade 平均 diff_axis:\n\n{by.to_markdown()}\n"]

    # ② stage 校正 Cox
    cd = d.dropna(subset=["time", "event", "diff_axis", "stage_ord", "age_hi"])
    cd = cd[cd["time"] > 0]
    cph = CoxPHFitter()
    cph.fit(cd[["time", "event", "diff_axis", "stage_ord", "age_hi"]], "time", "event")
    s = cph.summary
    L += ["## ② stage+age 校正多变量 Cox\n", f"- n={len(cd)}\n",
          s[["coef", "exp(coef)", "p"]].round(4).to_markdown() + "\n",
          f"- **diff_axis 校正后 HR={np.exp(s.loc['diff_axis','coef']):.3f}, p={s.loc['diff_axis','p']:.4f}**\n"]

    verdict = ("分化签名吻合病理分级（签名有效），但无独立预后价值"
               if pg < 0.05 and rho < 0 and s.loc["diff_axis", "p"] > 0.05
               else "见上方数字")
    L += [f"## 小结\n- {verdict}\n"]

    with open(f"{C.QC}/S6b_grade_adjust.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\ndiff_axis~Grade rho={rho:.3f} p={pg:.4f}; adj-Cox p={s.loc['diff_axis','p']:.4f}")


if __name__ == "__main__":
    main()
