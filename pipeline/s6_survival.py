#!/usr/bin/env python
"""s6_survival.py — 分化态的预后价值（bulk 生存分析）。★第 4 根支柱★

假设：低分化/高 M2（低 diff_axis）的肿瘤预后更差。在 TCGA-STAD 和 GSE84437 两套独立
bulk 验证 → 跨数据集稳健才算数。签名来自 S3_5_states.h5ad（与单细胞一致，去 MDK）。

依赖: pip install lifelines
Output: B8GC_WORK_ROOT/qc/S6_survival.md
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import config as C
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

META = f"{C.RESULTS}/S3_5_states.h5ad"
STATES = ["M0", "M1", "M2"]
MDK_AVOID = {"MDK"}

DATASETS = {
    "TCGA-STAD": dict(
        expr=f"{C.HUB}/GC/TCGA_STAD_Exp.txt", expr_sep=",",
        gene_col="gene_name", drop_cols=["gene_id", "gene_type"],
        tumor_suffix="-01", truncate_barcode=True,
        clin=f"{C.HUB}/GC/TCGA_STAD_time.txt", clin_sep="\t",
        sample_col="Id", time_col="futime", event_col="fustat",
    ),
    "GSE84437": dict(
        expr=f"{C.HUB}/GC/GSE84437_raw/GSE84437_estimate_input.txt", expr_sep="\t",
        gene_col="Gene", drop_cols=[],
        tumor_suffix=None, truncate_barcode=False,
        clin=f"{C.HUB}/GC/GSE84437_raw/GSE84437_survival_raw.csv", clin_sep=",",
        sample_col="sample", time_col="duration overall survival:ch1", event_col="death:ch1",
    ),
}


def derive_sigs():
    from sig_utils import frozen_sigs
    return frozen_sigs()


def load_expr(cfg):
    e = pd.read_csv(cfg["expr"], sep=cfg["expr_sep"])
    e = e.drop(columns=[c for c in cfg["drop_cols"] if c in e.columns], errors="ignore")
    e = e.set_index(cfg["gene_col"])
    e = e.apply(pd.to_numeric, errors="coerce")
    e = e.groupby(level=0).mean()                       # 重复 symbol 取均值
    if cfg["tumor_suffix"]:                              # 只留肿瘤样本
        e = e.loc[:, [c for c in e.columns if str(c).endswith(cfg["tumor_suffix"])]]
    if cfg["truncate_barcode"]:                          # TCGA-CG-4444-01 → TCGA-CG-4444
        e.columns = ["-".join(str(c).split("-")[:3]) for c in e.columns]
        e = e.loc[:, ~pd.Index(e.columns).duplicated()]
    if np.nanmax(e.values) > 100:                        # 线性 → log
        e = np.log1p(e)
    return e


def score_bulk(expr, genes):
    g = [x for x in genes if x in expr.index]
    if len(g) < 3:
        return None
    z = expr.loc[g]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0)


def main():
    sigs = derive_sigs()
    L = ["# S6 分化态预后价值（bulk 生存）\n",
         f"- M1 签名 {len(sigs.get('M1',[]))} 基因; M2 签名 {len(sigs.get('M2',[]))} 基因\n"]
    for name, cfg in DATASETS.items():
        if not (os.path.exists(cfg["expr"]) and os.path.exists(cfg["clin"])):
            L.append(f"## {name}: 文件缺失 → 跳过\n"); continue
        expr = load_expr(cfg)
        s_M1, s_M2 = score_bulk(expr, sigs.get("M1", [])), score_bulk(expr, sigs.get("M2", []))
        if s_M1 is None or s_M2 is None:
            L.append(f"## {name}: 签名覆盖不足 → 跳过\n"); continue
        df = pd.DataFrame({"diff_axis": s_M1 - s_M2, "s_M1": s_M1, "s_M2": s_M2})
        df.index = df.index.astype(str); df.index.name = "sample"; df = df.reset_index()

        clin = pd.read_csv(cfg["clin"], sep=cfg["clin_sep"])
        clin = clin.rename(columns={cfg["sample_col"]: "sample",
                                    cfg["time_col"]: "time", cfg["event_col"]: "event"})
        clin["sample"] = clin["sample"].astype(str)
        d = df.merge(clin[["sample", "time", "event"]], on="sample")
        d["time"] = pd.to_numeric(d["time"], errors="coerce")
        d["event"] = pd.to_numeric(d["event"], errors="coerce")
        d = d.dropna(subset=["time", "event", "diff_axis"])
        d = d[d["time"] > 0]
        if len(d) < 20:
            L.append(f"## {name}: 合并后仅 {len(d)} 例（匹配可能失败）→ 跳过\n"); continue

        cph = CoxPHFitter(); cph.fit(d[["time", "event", "diff_axis"]], "time", "event")
        hr = float(np.exp(cph.params_["diff_axis"])); p = float(cph.summary.loc["diff_axis", "p"])
        d["grp"] = np.where(d["diff_axis"] > d["diff_axis"].median(), "high_diff", "low_diff")
        lr = logrank_test(d[d.grp == "high_diff"]["time"], d[d.grp == "low_diff"]["time"],
                          d[d.grp == "high_diff"]["event"], d[d.grp == "low_diff"]["event"])
        verdict = ("低分化(高M2)预后差" if hr < 1 and p < 0.05
                   else "高分化预后差" if hr > 1 and p < 0.05 else "无显著关联")
        L += [f"## {name}（n={len(d)}）\n",
              f"- Cox(diff_axis 连续): HR={hr:.3f} (>1=高分化差; <1=低分化差), p={p:.4f}",
              f"- KM(diff_axis 中位分组): log-rank p={lr.p_value:.4f}",
              f"- **解读: {verdict}**\n"]
        print(f"{name}: n={len(d)} HR={hr:.3f} cox_p={p:.4f} km_p={lr.p_value:.4f} -> {verdict}")

    with open(f"{C.QC}/S6_survival.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
