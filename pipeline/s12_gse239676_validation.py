#!/usr/bin/env python
"""s12_gse239676_validation.py — frozen-signature validation on GSE239676.

Uses the existing B8 frozen M0/M1/M2 signatures without redefining states.
Reports route-resolved diff_axis and intrinsic-8 EMT effects, plus optional
CopyKAT malignant-only sensitivity if sample sizes allow.
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp
from scipy.stats import pearsonr, mannwhitneyu, wilcoxon, zscore
import config as C
from s12_external_validation import epithelial, score_modules, run_copykat_on_epithelial, cohend
from sig_utils import frozen_sigs, signature_sha256, stamp_signature_hash

GSE239676_ROOT = os.environ.get("B8GC_GSE239676_ROOT", f"{C.HUB}/GSE239676_raw")
MTX = f"{GSE239676_ROOT}/GSE239676_count_matrix.mtx.gz"
BAR = f"{GSE239676_ROOT}/GSE239676_barcodes.tsv.gz"
FEAT = f"{GSE239676_ROOT}/GSE239676_features.tsv.gz"
META = f"{GSE239676_ROOT}/GSE239676_meta.tsv.gz"
COPYKAT_ROOT = os.environ.get("B8GC_COPYKAT_SCRATCH_239676", f"{C.WORK}/tmp/copykat_gse239676")
COPYKAT_MIN_CELLS = 100

TISSUE_TO_ROUTE = {
    "P": "primary",
    "Li": "liver_met",
    "As": "ascites",
    "Ov": "ovarian_met",
    "Ad": "adjacent_normal",
    "PB": "drop",
    "PBMC": "drop",
}
MET_ROUTES = ["liver_met", "ascites", "ovarian_met"]


def is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def read_table_auto(path, sep="\t"):
    return pd.read_csv(path, sep=sep, compression="gzip" if is_gzip(path) else None)


def read_lines_auto(path):
    if is_gzip(path):
        import gzip
        with gzip.open(path, "rt") as fh:
            return [x.rstrip("\n") for x in fh]
    with open(path, "rt") as fh:
        return [x.rstrip("\n") for x in fh]


def load_mtx_with_meta():
    p = MTX
    if is_gzip(MTX):
        import gzip, tempfile, shutil
        t = tempfile.NamedTemporaryFile(delete=False, suffix=".mtx")
        with gzip.open(MTX, "rb") as fi:
            shutil.copyfileobj(fi, t)
        t.close()
        p = t.name
    X = sio.mmread(p).T.tocsr()
    if p != MTX:
        os.unlink(p)
    barcodes = read_lines_auto(BAR)
    feats = [x.split("\t") for x in read_lines_auto(FEAT)]
    genes = [(r[1] if len(r) > 1 else r[0]) for r in feats]
    meta = read_table_auto(META, sep="\t")
    if len(barcodes) != meta.shape[0]:
        raise ValueError(f"barcode rows {len(barcodes)} != meta rows {meta.shape[0]}")
    ad = sc.AnnData(X)
    ad.obs_names = pd.Index(barcodes[:ad.n_obs])
    ad.var_names = pd.Index(genes[:ad.n_vars])
    ad.var_names_make_unique()
    ad.obs = meta.iloc[:ad.n_obs].copy().reset_index(drop=True)
    ad.obs.index = ad.obs_names
    ad.obs["Sample"] = ad.obs["Sample"].astype(str)
    ad.obs["Patient"] = ad.obs["Patient"].astype(str)
    ad.obs["Tissue"] = ad.obs["Tissue"].astype(str)
    ad.obs["sample"] = ad.obs["Sample"]
    ad.obs["patient"] = ad.obs["Patient"]
    ad.obs["route"] = ad.obs["Tissue"].map(TISSUE_TO_ROUTE).fillna("unknown")
    ad.obs["cohort"] = "GSE239676"
    ad.obs["barcode"] = ad.obs_names.astype(str)
    return ad


def sample_level(df, route):
    pri = (df[df["route"] == "primary"]
             .groupby("sample", observed=True)
             .agg(dax=("dax", "mean"), emt=("s_EMT", "mean"), emti=("s_EMTi", "mean"), patient=("patient", "first"))
             .reset_index())
    met = (df[df["route"] == route]
             .groupby("sample", observed=True)
             .agg(dax=("dax", "mean"), emt=("s_EMT", "mean"), emti=("s_EMTi", "mean"), patient=("patient", "first"))
             .reset_index())
    return pri, met


def paired_delta(pri, met, value_col):
    p = pri[["patient", value_col]].rename(columns={value_col: "primary"})
    m = met[["patient", value_col]].rename(columns={value_col: "met"})
    sub = p.merge(m, on="patient", how="inner")
    if sub.shape[0] == 0:
        return sub, np.nan, np.nan
    delta = sub["met"] - sub["primary"]
    pval = np.nan
    if sub.shape[0] >= 2 and not np.allclose(sub["primary"], sub["met"]):
        try:
            pval = wilcoxon(sub["primary"], sub["met"]).pvalue
        except Exception:
            pval = np.nan
    return sub.assign(delta=delta), float(delta.mean()), pval


def route_report(df, route, L):
    pri = df.loc[df.route == "primary", "dax"]
    met = df.loc[df.route == route, "dax"]
    if len(pri) < 20 or len(met) < 20:
        L.append(f"## {route}\n- 细胞不足: primary={len(pri)}, {route}={len(met)}\n")
        return
    delta = float(met.mean() - pri.mean())
    pri_s, met_s = sample_level(df, route)
    L.append(f"## {route}")
    L.append(f"- 细胞级 diff_axis: primary {pri.mean():+.3f} → {route} {met.mean():+.3f}; Δ={delta:+.3f}")
    L.append(f"- 样本级均值: primary={list(pri_s['dax'].round(2))} vs {route}={list(met_s['dax'].round(2))}")
    try:
        _, ps = mannwhitneyu(pri_s["dax"], met_s["dax"], alternative="greater")
    except Exception:
        ps = np.nan
    try:
        _, pc = mannwhitneyu(pri, met, alternative="greater")
    except Exception:
        pc = np.nan
    L.append(f"- Mann-Whitney: sample-level p={ps:.3g} (n={pri_s.shape[0]}v{met_s.shape[0]}), cell-level p={pc:.3g} (描述性)")
    pair_dax, mean_dd, p_pair = paired_delta(pri_s, met_s, "dax")
    if pair_dax.shape[0] > 0:
        L.append(f"- 配对病人 Δdiff_axis: n={pair_dax.shape[0]}, mean Δ={mean_dd:+.3f}, Wilcoxon p={p_pair:.3g}")
    dd_diff = cohend(pri, met)
    ep = df.loc[df.route == "primary", "s_EMT"]
    em = df.loc[df.route == route, "s_EMT"]
    epi = df.loc[df.route == "primary", "s_EMTi"]
    emi = df.loc[df.route == route, "s_EMTi"]
    dd_emt = cohend(ep, em)
    dd_emti = cohend(epi, emi)
    L.append(f"- Cohen's d: diff_axis={dd_diff:+.2f}, EMT_full={dd_emt:+.2f}, EMT_intrinsic8={dd_emti:+.2f}")
    p_emt, mean_emt, p_emt_pair = paired_delta(pri_s, met_s, "emti")
    if p_emt.shape[0] > 0:
        L.append(f"- 配对病人 ΔEMT_intrinsic8: n={p_emt.shape[0]}, mean Δ={mean_emt:+.3f}, Wilcoxon p={p_emt_pair:.3g}")
    if abs(dd_emti) >= abs(dd_diff):
        verdict = "⚠ intrinsic-8 EMT 至少与去分化同量级，不能写 non-EMT"
    elif abs(dd_emti) >= 0.4:
        verdict = "去分化为主，但有真实 intrinsic-8 EMT 成分"
    else:
        verdict = "intrinsic-8 EMT 较小，去分化主导"
    L.append(f"- 判读: {verdict}\n")


def lauren_report(df, meta_cols, L):
    cands = [c for c in meta_cols if any(k in c.lower() for k in ["lauren", "hist", "diffuse", "intestinal"])]
    if not cands:
        L.append("## Lauren\n- metadata 无 Lauren / histology 列，无法评估。\n")
        return
    L.append(f"## Lauren\n- 候选列: {cands}\n")


def main(copykat=False, copykat_cores=8):
    sigs = frozen_sigs()
    sig_hash = signature_sha256(sigs)
    L = ["# S12 GSE239676 frozen-signature validation"]
    ad = load_mtx_with_meta()
    L.append(f"- 原始载入: {ad.n_obs} cells × {ad.n_vars} genes")
    L.append(f"- metadata 列: {list(ad.obs.columns)}")
    L.append(f"- Tissue 分布: {dict(ad.obs['Tissue'].value_counts())}")
    tissue_map = (ad.obs[["sample", "patient", "Tissue", "route"]]
                    .drop_duplicates()
                    .sort_values(["patient", "sample"]))
    L.append("- sample × Tissue × route 映射:\n" + tissue_map.to_string(index=False))
    ad = ad[~ad.obs["route"].isin(["drop", "unknown"])].copy()
    L.append(f"- 去掉 blood/unknown 后: {ad.n_obs} cells; route 分布 {dict(ad.obs['route'].value_counts())}")

    e, cl, keep = epithelial(ad)
    L.append(f"- 上皮簇 {keep} → epithelial {e.n_obs} cells")
    L.append(f"- 上皮簇分数: {cl.round(3).to_dict()}")

    if copykat:
        cfg = {"name": "GSE239676"}
        emal, manifest = run_copykat_on_epithelial(e, cfg, COPYKAT_ROOT, ncores=copykat_cores, min_cells=COPYKAT_MIN_CELLS)
        manifest.to_csv(f"{C.RESULTS}/S12_239676_copykat_manifest.csv", index=False)
        L.append(f"- CopyKAT manifest: {C.RESULTS}/S12_239676_copykat_manifest.csv")
        for _, rr in manifest.iterrows():
            L.append(f"  · {rr['sample']} [{rr['route']}]: epi={int(rr['n_epi'])} pred={int(rr['n_pred'])} aneu={int(rr['n_aneuploid'])} status={rr['status']}")
        if emal.n_obs > 0:
            e = emal
            L.append(f"- CopyKAT malignant-only retained {e.n_obs} cells")
        else:
            L.append("- ⚠ CopyKAT 后无 malignant 细胞，回退 epithelial-gated 主分析")

    emt_src, nEMT, nEMTtot, nEMTi = score_modules(e, sigs)
    L.append(f"- M1 hit {len([x for x in sigs['M1'] if x in e.var_names])}/{len(sigs['M1'])}; M2 hit {len([x for x in sigs['M2'] if x in e.var_names])}/{len(sigs['M2'])}")
    L.append(f"- EMT source: {emt_src}; full {nEMT}/{nEMTtot}; intrinsic-8 {nEMTi}")
    r, p = pearsonr(e.obs["s_M1"], e.obs["s_M2"])
    L.append(f"- Axis reproduced: M1 vs M2 r={r:.3f}, p={p:.2e}")

    out = e.obs[["sample", "patient", "route", "s_M1", "s_M2", "dax", "s_EMT", "s_EMTi"]].copy()
    out["cohort"] = "GSE239676"
    out = stamp_signature_hash(out, sig_hash)
    out.to_csv(f"{C.RESULTS}/S12_239676_cells.csv", index=False)
    L.append(f"- per-cell output: {C.RESULTS}/S12_239676_cells.csv")

    for route in MET_ROUTES:
        if (e.obs["route"] == route).sum() > 0:
            route_report(out, route, L)

    lauren_report(out, list(ad.obs.columns), L)
    with open(f"{C.QC}/S12_239676_external.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--copykat", action="store_true", help="Run sample-level CopyKAT sensitivity")
    ap.add_argument("--copykat-cores", type=int, default=8)
    args = ap.parse_args()
    main(copykat=args.copykat, copykat_cores=args.copykat_cores)
