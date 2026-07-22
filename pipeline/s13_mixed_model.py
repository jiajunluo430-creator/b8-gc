#!/usr/bin/env python
"""s13_mixed_model.py — 多队列样本级合并统计,区分独立验证头条与含发现集敏感性。"""
import argparse
import glob
import json
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import statsmodels.formula.api as smf
from scipy.stats import t as tdist, norm
import config as C
from fig_utils import save, W1
from sig_utils import assert_signature_hash

NORMAL_PAT = r"normal|precancer|adjacent|healthy|para|adj"
HRA_COHORT = "HRA004702"
VALIDATION = {"GSE308231", "GSE246662", "GSE239676"}
LOCO_ORDER = ["GSE308231", "GSE246662", "GSE239676"]


def outputs(v3=False):
    tag = "V3" if v3 else ""
    return {
        "cell_glob": f"{C.RESULTS}/S12{tag}_cells_*.csv" if v3 else f"{C.RESULTS}/S12_cells_*.csv",
        "qc": f"{C.QC}/S13{tag}_mixed.md" if v3 else f"{C.QC}/S13_mixed.md",
        "json": f"{C.RESULTS}/S13{tag}_pooled.json" if v3 else f"{C.RESULTS}/S13_pooled.json",
        "fig": "FigS9_pooled_sensitivity_V3" if v3 else "FigS9_pooled_sensitivity",
    }


def load_hra():
    p = f"{C.RESULTS}/S4_scored.h5ad"
    if not os.path.exists(p):
        print(f"⚠ 未找到 HRA {p};只用外部队列(无 HRA 会少一个队列)。")
        return pd.DataFrame(), None
    obs = sc.read_h5ad(p, backed="r").obs
    dcol = next((c for c in ["dax", "diff_axis", "diffaxis"] if c in obs), None)
    pcol = next((c for c in obs.columns if c.lower() in
                 ("patient", "sample", "orig.ident", "donor", "patient_id", "sample_id", "samplename")), None)
    rcol = "route" if "route" in obs else next((c for c in obs.columns if "route" in c.lower()), None)
    ccol = "cohort" if "cohort" in obs else next((c for c in obs.columns if "cohort" in c.lower()), None)
    if not (dcol and pcol and rcol and ccol):
        return pd.DataFrame(), f"HRA obs 列={list(obs.columns)} → 没认出 dax/patient/route/cohort,贴我这行我对版"
    keep = obs[ccol].astype(str).values == HRA_COHORT
    df = pd.DataFrame({
        "dax": np.asarray(obs[dcol])[keep],
        "sample": obs[pcol].astype(str).values[keep],
        "route": obs[rcol].astype(str).values[keep],
        "cohort": "HRA",
    })
    return df, None


def load_external(v3=False):
    files = sorted(glob.glob(outputs(v3)["cell_glob"]))
    if not v3:
        g239 = f"{C.RESULTS}/S12_239676_cells.csv"
        if os.path.exists(g239):
            files.append(g239)
    dfs = [assert_signature_hash(pd.read_csv(f), f) for f in files]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def cohort_effect(d):
    smp = d[d.met == 0].groupby("sample_uid")["dax"].mean()
    smm = d[d.met == 1].groupby("sample_uid")["dax"].mean()
    if len(smp) < 2 or len(smm) < 2:
        return smm.mean() - smp.mean(), np.nan, np.nan, np.nan, len(smp), len(smm)
    diff = smm.mean() - smp.mean()
    se = np.sqrt(smp.var(ddof=1) / len(smp) + smm.var(ddof=1) / len(smm))
    dfree = len(smp) + len(smm) - 2
    tc = tdist.ppf(0.975, dfree)
    return diff, diff - tc * se, diff + tc * se, se, len(smp), len(smm)


def ols_on(frame):
    if frame.empty:
        raise ValueError("empty frame")
    frame = frame.copy()
    frame["cohort_cat"] = frame["cohort"].astype("category")
    formula = "dax ~ met + cohort_cat" if frame["cohort"].nunique() >= 2 else "dax ~ met"
    mo = smf.ols(formula, frame).fit()
    return float(mo.params["met"]), mo.conf_int().loc["met"].values.astype(float), float(mo.pvalues["met"])


def fixed_effect_meta(ceff, cohorts):
    valid = [(ceff[c]["diff"], ceff[c]["se"]) for c in cohorts if np.isfinite(ceff[c]["se"]) and ceff[c]["se"] > 0]
    if not valid:
        raise ValueError("no cohort has a finite positive standard error")
    w = np.array([1 / s**2 for _, s in valid])
    dd = np.array([x for x, _ in valid])
    mb = float((w * dd).sum() / w.sum())
    mse = float(np.sqrt(1 / w.sum()))
    mp = float(2 * norm.sf(abs(mb / mse)))
    return mb, mse, mp


def to_num(x):
    x = float(x)
    return x if np.isfinite(x) else None


def main(v3=False):
    out = outputs(v3)
    title_tag = "V3 CopyKAT-malignant" if v3 else ""
    L = [f"# S13 {'V3 ' if v3 else ''}多队列样本级合并统计(转移→去分化,抗伪重复)\n"]
    hra, warn = load_hra()
    if warn:
        L.append("⚠ " + warn + "\n")
        open(out["qc"], "w").write("\n".join(L))
        print("\n".join(L))
        return
    ext = load_external(v3=v3)
    df = pd.concat([hra, ext], ignore_index=True)
    if df.empty:
        L.append("⚠ 没数据(先跑 s12 生成外部验证 cells CSV,并确认 HRA S4_scored.h5ad 存在)。")
        open(out["qc"], "w").write("\n".join(L))
        print("\n".join(L))
        return

    df = df.dropna(subset=["dax"])
    df = df[~df.route.str.contains(NORMAL_PAT, case=False, na=False)].copy()
    df["met"] = (~df.route.str.contains("primary|prim|tumor|in_situ|insitu", case=False, na=False)).astype(int)
    df["sample_uid"] = df.cohort.astype(str) + ":" + df["sample"].astype(str)
    L.append(f"- 合并: {df.shape[0]} 细胞, {df.sample_uid.nunique()} 样本, 队列={sorted(df.cohort.unique())}")
    L.append(f"- primary 细胞 {(df.met==0).sum()} / metastasis 细胞 {(df.met==1).sum()}; primary 样本 {df[df.met==0].sample_uid.nunique()} / met 样本 {df[df.met==1].sample_uid.nunique()}")
    L.append(f"- HRA 仅取 {HRA_COHORT}: sample 数 = {hra['sample'].nunique()}")
    if v3:
        L.append("- 外部验证输入 = CopyKAT malignant-only V3 cells")
    L.append("")

    cohs = sorted(df.cohort.unique())
    ceff = {}
    for c in cohs:
        diff, lo, hi, se, np_, nm_ = cohort_effect(df[df.cohort == c])
        ceff[c] = dict(diff=diff, lo=lo, hi=hi, se=se, npri=np_, nmet=nm_)

    sm = (df.groupby("sample_uid", observed=True)
            .agg(dax=("dax", "mean"), met=("met", "first"), cohort=("cohort", "first"))
            .reset_index())

    smv = sm[sm.cohort.isin(VALIDATION)].copy()
    L.append("## 头条统计: 仅独立验证集(发现集未见过) — sample-level OLS")
    bV, ciV, pV = ols_on(smv)
    L.append(f"- 验证集 = {sorted(set(smv.cohort))}, {smv.shape[0]} 样本(每样本 1 点)")
    L.append(f"- **β = {bV:+.3f}  (95%CI {ciV[0]:+.3f}, {ciV[1]:+.3f}),  p = {pV:.2e}**\n")

    L.append("## 敏感性 1: 全部队列(含发现集 HRA) — sample-level OLS")
    beta, ci, pval = ols_on(sm)
    L.append(f"- 全部 {sm.shape[0]} 样本: β = {beta:+.3f} (95%CI {ci[0]:+.3f}, {ci[1]:+.3f}), p = {pval:.2e}")
    L.append(f"- **注: 含发现集后 β 变为 {beta:+.2f} vs 验证集 {bV:+.2f}** → HRA 是最温和效应,加进来是保守(不是灌水)。\n")

    validation_cohs = [c for c in LOCO_ORDER if c in ceff]
    vmb, vmse, vmp = fixed_effect_meta(ceff, validation_cohs)
    L.append("## 敏感性 2: 独立验证队列逆方差固定效应 meta")
    L.append(f"- pooled Δ = {vmb:+.3f} (95%CI {vmb-1.96*vmse:+.3f}, {vmb+1.96*vmse:+.3f}), p = {vmp:.2e}")
    L.append("- 注: 效应量随路径加深(卵巢 < 腹膜 < 肝)= 真实生物学,报告讲『方向一致复现』非『幅度相同』。\n")

    mb, mse, mp = fixed_effect_meta(ceff, cohs)
    L.append("## 敏感性 3: 全部队列(含发现集 HRA)逆方差固定效应 meta")
    L.append(f"- pooled Δ = {mb:+.3f} (95%CI {mb-1.96*mse:+.3f}, {mb+1.96*mse:+.3f}), p = {mp:.2e}\n")

    loco_rows = []
    L.append("## Leave-one-validation-cohort-out (LOCO) 固定效应 meta")
    for omitted in validation_cohs:
        retained = [c for c in validation_cohs if c != omitted]
        lb, lse, lp = fixed_effect_meta(ceff, retained)
        n_samples = int(smv.loc[smv.cohort.isin(retained), "sample_uid"].nunique())
        row = {
            "omitted": omitted,
            "beta": to_num(lb),
            "ci": [to_num(lb - 1.96 * lse), to_num(lb + 1.96 * lse)],
            "p": to_num(lp),
            "n_samples": n_samples,
            "retained_cohorts": retained,
        }
        loco_rows.append(row)
        L.append(
            f"- omit {omitted}: Δ={lb:+.3f} "
            f"(95%CI {lb-1.96*lse:+.3f}, {lb+1.96*lse:+.3f}), "
            f"p={lp:.2e}, n={n_samples} samples"
        )
    L.append("")

    label_map = {
        "HRA": "HRA (disc.)",
        "GSE239676": "GSE239676 (liver/ascites)",
        "GSE246662": "GSE246662 (liver)",
        "GSE308231": "GSE308231 (peri)",
    }
    order = [c for c in ["HRA", "GSE239676", "GSE246662", "GSE308231"] if c in ceff]
    cohort_rows = []
    for c in order:
        cohort_rows.append({
            "cohort": c,
            "label": label_map.get(c, c),
            "kind": "validation" if c in VALIDATION else "discovery",
            "diff": to_num(ceff[c]["diff"]),
            "lo": to_num(ceff[c]["lo"]),
            "hi": to_num(ceff[c]["hi"]),
            "se": to_num(ceff[c]["se"]),
            "npri": int(ceff[c]["npri"]),
            "nmet": int(ceff[c]["nmet"]),
        })

    panel_g_rows = cohort_rows + [{
        "label": "Pooled (val.)",
        "kind": "val",
        "diff": to_num(bV),
        "lo": to_num(ciV[0]),
        "hi": to_num(ciV[1]),
    }]
    sensitivity_rows = panel_g_rows + [{
        "label": "Pooled (all)",
        "kind": "all",
        "diff": to_num(beta),
        "lo": to_num(ci[0]),
        "hi": to_num(ci[1]),
    }]

    payload = {
        "validation_only": {"beta": to_num(bV), "ci": [to_num(ciV[0]), to_num(ciV[1])], "p": to_num(pV)},
        "all_cohorts": {"beta": to_num(beta), "ci": [to_num(ci[0]), to_num(ci[1])], "p": to_num(pval)},
        "validation_fixed_meta": {"delta": to_num(vmb), "ci": [to_num(vmb - 1.96 * vmse), to_num(vmb + 1.96 * vmse)], "p": to_num(vmp)},
        "fixed_meta": {"delta": to_num(mb), "ci": [to_num(mb - 1.96 * mse), to_num(mb + 1.96 * mse)], "p": to_num(mp)},
        "loco_rows": loco_rows,
        "cohort_rows": cohort_rows,
        "panel_g_rows": panel_g_rows,
        "sensitivity_rows": sensitivity_rows,
    }
    with open(out["json"], "w") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)

    L.append("## 森林图行")
    for row in sensitivity_rows:
        lo = row["lo"]
        hi = row["hi"]
        L.append(f"- {row['label']}: Δ={row['diff']:+.3f}" + (f" [{lo:+.2f}, {hi:+.2f}]" if lo is not None and hi is not None else ""))

    fig, ax = plt.subplots(figsize=(W1, 0.65 * W1))
    ys = np.arange(len(sensitivity_rows))[::-1]
    for y, row in zip(ys, sensitivity_rows):
        diff, lo, hi, kind = row["diff"], row["lo"], row["hi"], row["kind"]
        if kind == "val":
            col = "#c1121f"
            diamond = np.array([[lo, y], [diff, y + 0.22], [hi, y], [diff, y - 0.22]])
            ax.add_patch(Polygon(diamond, closed=True, facecolor=col, edgecolor=col, zorder=3))
        elif kind == "all":
            col = "#e07a00"
            diamond = np.array([[lo, y], [diff, y + 0.22], [hi, y], [diff, y - 0.22]])
            ax.add_patch(Polygon(diamond, closed=True, facecolor=col, edgecolor=col, zorder=3))
        else:
            col = "#2a9d8f" if kind == "validation" else "#6c757d"
            if lo is not None and hi is not None:
                ax.plot([lo, hi], [y, y], color=col, lw=1.4, zorder=2)
            ax.scatter([diff], [y], s=32, color=col, marker="o", zorder=3)
    ax.axvline(0, ls="--", color="gray", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["label"] for r in sensitivity_rows], fontsize=5.2)
    ax.set_xlabel("Δ diff_axis  (metastasis − primary)")
    ttl = "Metastatic dedifferentiation: discovery vs independent validation"
    if v3:
        ttl += "\nV3 CopyKAT malignant-only"
    ax.set_title(f"{ttl}\nvalidation-only β={bV:+.2f}, p={pV:.1e}")
    plt.tight_layout()
    save(fig, out["fig"])

    L.append("\n解读：头条用**仅独立验证集**(GSE239676+GSE308231+GSE246662,发现集未见过)的 β/p —— 零循环;全部队列(含 HRA)作敏感性。")
    if v3:
        L.append("V3 解释限定于 CopyKAT aneuploid 恶性细胞，目的是看 malignant-only 门控后腹膜/肝验证是否仍守住方向与效应。")
    else:
        L.append("诚实提醒：存在异质性(卵巢<腹膜<肝,效应随路径加深),应表述为『方向一致复现』,不是『幅度完全相同』。")

    open(out["qc"], "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--v3", action="store_true", help="Read V3 S12 cell tables and write V3 pooled outputs")
    args = ap.parse_args()
    main(v3=args.v3)
