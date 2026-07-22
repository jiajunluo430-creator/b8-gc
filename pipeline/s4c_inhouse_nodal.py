#!/usr/bin/env python
"""s4c_inhouse_nodal.py — 自采原发灶 N+ vs N0 去分化(临床关联:去分化预示淋巴结转移倾向)。

自采 3v3 = 6 个原发灶，3 个淋巴结阳性(N+)、3 个阴性(N0)。
测:原发灶里 diff_axis / M2 占比 是否在 N+ 更去分化(与"去分化预示转移倾向"一致)。
n=3v3 极小 → 方向性/探索性证据，非主统计。

先自检淋巴结状态如何编码(obs 列 / 样本名)，能识别就出对比，识别不了就打印诊断让我对版。
Output: B8GC_WORK_ROOT/qc/S4c_inhouse_nodal.md
"""
import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
import config as C

SC = f"{C.RESULTS}/S4_scored.h5ad"


def detect_nodal(ih):
    """尝试从 obs 列或样本名识别 N+/N0。返回 Series(index=cell, 值∈{'N+','N0'}) 或 None。"""
    # 1) 直接有淋巴结状态列
    for c in ih.columns:
        vals = ih[c].astype(str).str.upper()
        uset = set(vals.unique())
        if uset & {"N+", "N1", "N2", "N3", "POSITIVE", "POS"} and uset & {"N0", "NEGATIVE", "NEG"}:
            return vals.map(lambda v: "N+" if v in {"N+", "N1", "N2", "N3", "POSITIVE", "POS"} else "N0"), f"列 {c}"
        if c.lower() in {"ln_status", "node", "nodal", "n_stage", "ln", "lnm", "metastasis"}:
            return vals.map(lambda v: "N0" if v in {"0", "N0", "NEG", "NEGATIVE", "FALSE", "NO"} else "N+"), f"列 {c}"
    # 2) 自采 3v3 固定样本映射
    s = ih["sample"].astype(str)
    inhouse_map = {
        "P01T": "N+", "P02T": "N+", "P04T": "N+",
        "P07T": "N0", "P12T": "N0", "P19T": "N0",
    }
    if set(s.unique()).issubset(set(inhouse_map)):
        return s.map(inhouse_map), "自采样本映射(P01/P02/P04=N+; P07/P12/P19=N0)"
    # 3) 样本名里含 N0/N1 或 LN/pos/neg
    if s.str.upper().str.contains(r"N[+1-3]").any() and s.str.upper().str.contains("N0").any():
        return s.map(lambda x: "N0" if re.search(r"N0", x.upper()) else ("N+" if re.search(r"N[+1-3]", x.upper()) else np.nan)), "样本名 N0/N1"
    if s.str.lower().str.contains("pos|lnp|n_pos").any() and s.str.lower().str.contains("neg|lnn|n_neg").any():
        return s.map(lambda x: "N+" if re.search("pos|lnp", x.lower()) else "N0"), "样本名 pos/neg"
    return None, None


def main():
    a = sc.read_h5ad(SC); o = a.obs.copy()
    if "meta" not in o and "state" in o:
        o["meta"] = o["state"]
    o["meta"] = o["meta"].astype(str)
    o["cohort"] = o["cohort"].astype(str)
    ih = o[o["cohort"].str.upper().str.contains("INHOUSE|3V3")].copy()
    L = ["# S4c 自采原发灶 N+ vs N0 去分化(临床关联)\n"]
    if ih.empty:
        L.append("⚠ 未找到 INHOUSE/3v3 队列；贴回 S4_scored 的 cohort 取值。")
        open(f"{C.QC}/S4c_inhouse_nodal.md", "w").write("\n".join(L)); print("\n".join(L)); return

    # 诊断信息(识别不了时给我对版用)
    L += ["## 自检", f"- 自采细胞 {len(ih)}；obs 列={list(ih.columns)}",
          f"- sample 取值={sorted(ih['sample'].astype(str).unique())}",
          f"- route 取值={dict(ih['route'].value_counts())}\n"]

    grp, src = detect_nodal(ih)
    if grp is None or grp.isna().all():
        L.append("⚠ **未能自动识别 N+/N0**。请告诉我淋巴结状态在哪个列/或样本名规则，我对版改。")
        open(f"{C.QC}/S4c_inhouse_nodal.md", "w").write("\n".join(L)); print("\n".join(L)); return
    ih["nodal"] = grp.values
    L.append(f"- ✓ 识别 N+/N0 来源：{src}；分组 n={dict(ih['nodal'].value_counts())}\n")

    npos = ih[ih.nodal == "N+"]; nneg = ih[ih.nodal == "N0"]
    # per-cell 描述
    L += ["## ① per-cell diff_axis（描述性）",
          f"- N0 均值 {nneg['diff_axis'].mean():.3f} (n={len(nneg)}) → N+ 均值 {npos['diff_axis'].mean():.3f} (n={len(npos)})",
          f"- Δ(N+ − N0)={npos['diff_axis'].mean()-nneg['diff_axis'].mean():+.3f}（负=N+ 更去分化，符合预期）\n"]
    # per-sample pseudobulk（主依据）
    psb = ih.groupby(["sample", "nodal"], observed=True)["diff_axis"].mean().reset_index()
    pp = psb[psb.nodal == "N+"]["diff_axis"]; pn = psb[psb.nodal == "N0"]["diff_axis"]
    L += ["## ② per-sample pseudobulk（主依据；n=3v3，看方向）",
          f"- N+ 样本 n={len(pp)} 均值 {pp.mean():.3f}；N0 样本 n={len(pn)} 均值 {pn.mean():.3f}",
          f"- Δ={pp.mean()-pn.mean():+.3f}"]
    if len(pp) >= 2 and len(pn) >= 2:
        u, pv = mannwhitneyu(pn, pp, alternative="greater")   # H1: N0 > N+ (N+ 更去分化)
        L.append(f"- Mann-Whitney(N0>N+) p={pv:.3f}（n=3v3 最小可达 p≈0.05–0.1，方向为主）\n")
    # M2 占比
    m2 = ih.groupby("nodal", observed=True)["meta"].apply(lambda s: (s == "M2").mean())
    L += ["## ③ M2(去分化)占比", m2.round(3).to_markdown(), "\n",
          "解读：若 N+ 原发灶 diff_axis 更低、M2 占比更高 → **原发灶去分化预示淋巴结转移倾向**，"
          "与 HRA 的'转移伴随去分化'互补，给纯干稿补上一条自采数据的临床关联。n 小，按探索性写。"]
    open(f"{C.QC}/S4c_inhouse_nodal.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
