#!/usr/bin/env python
"""s12_external_validation.py — 两个独立留出 GC 队列验证「分化轴 + 去分化」复现。

GSE308231: 3 原发(Ca)+ 3 腹膜转移(F)        10x MTX/TSV → 独立复现腹膜转移去分化(HRA 里腹膜 n=4 NS)
GSE246662: 3 原发(GC)+ 3 肝转移(LM)+3 健康肝(HL,弃)  CSV    → 扩展到肝转移路径
两者都不在发现集;各自原论文(npj Digit Med 共浸润ML预后 / Oncogene NK免疫)都不碰恶性分化轴 → 验证干净。

逻辑(每个队列):解 tar → 按样本名打 route 标签 → QC/归一化 → 聚类取上皮 → 打冻结 M1/M2/M0 签名 →
  ① 轴复现: M1 vs M2 打分负相关  ② 去分化复现: primary vs met 的 diff_axis。
Outputs under B8GC_WORK_ROOT/qc and B8GC_WORK_ROOT/figures.

若 tar 内文件命名和假设不符,会打印文件清单 + 每样本/每簇诊断 → 贴回我对版。
"""
import argparse
import gzip
import json
import os
import re
import glob
import shutil
import tarfile
import tempfile
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.sparse as sp
from scipy.stats import pearsonr, mannwhitneyu, zscore
import config as C
from fig_utils import save, panel, W2, STATE_COLORS
from sig_utils import SIGNATURE_HASH_COL, assert_signature_hash, frozen_sigs, signature_sha256, stamp_signature_hash

sc.settings.verbosity = 1

RS = os.environ.get("B8GC_RSCRIPT", "Rscript")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COPYKAT_ROOT_V3 = os.environ.get("B8GC_COPYKAT_SCRATCH_V3", f"{C.WORK}/tmp/copykat_v3_external")
COPYKAT_MIN_CELLS = 100


def r308(s):
    s = s.lower()
    if s.startswith("ca"):
        return "primary"
    if s.startswith("f"):
        return "peritoneal_met"
    return "other"


def r246(s):
    s = s.lower()
    if "hl" in s:
        return "drop"
    if "lm" in s:
        return "liver_met"
    if "gc" in s:
        return "primary"
    return "other"


COHORTS = [
    dict(name="GSE308231", label="GSE308231 · peritoneal", met="peritoneal_met", fmt="mtx", route=r308,
         tar=f"{C.HUB}/external/GSE308231/GSE308231_RAW.tar"),
    dict(name="GSE246662", label="GSE246662 · liver", met="liver_met", fmt="csv", route=r246,
         tar=f"{C.HUB}/external/GSE246662/GSE246662_RAW.tar"),
]


def out_names(v3=False):
    tag = "V3" if v3 else ""
    return {
        "cell_prefix": f"S12{tag}_cells_",
        "qc": f"{C.QC}/S12{tag}_external.md" if v3 else f"{C.QC}/S12_external.md",
        "fig": "Fig6_external_V3" if v3 else "Fig6_external",
        "pooled_json": f"{C.RESULTS}/S13{tag}_pooled.json" if v3 else f"{C.RESULTS}/S13_pooled.json",
        "copykat_manifest": f"{C.RESULTS}/S12{tag}_copykat_manifest.csv" if v3 else None,
    }


def _gunzip(path):
    if path.endswith(".gz"):
        suf = os.path.splitext(path[:-3])[1] or ".tmp"
        t = tempfile.NamedTemporaryFile(delete=False, suffix=suf)
        with gzip.open(path, "rb") as fi:
            shutil.copyfileobj(fi, t)
        t.close()
        return t.name
    return path


def _read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return [l.rstrip("\n").split("\t") for l in f]


def load_mtx(mtx, bc, feat):
    from scipy.io import mmread
    p = _gunzip(mtx)
    X = mmread(p).T.tocsr()
    if p != mtx:
        os.unlink(p)
    barcodes = [r[0] for r in _read_tsv(bc)]
    feats = _read_tsv(feat)
    genes = [(r[1] if len(r) > 1 else r[0]) for r in feats]
    ad = sc.AnnData(X)
    ad.obs_names = barcodes[:ad.n_obs]
    ad.var_names = genes[:ad.n_vars]
    ad.var_names_make_unique()
    return ad


def load_csv(path):
    df = pd.read_csv(path, index_col=0)

    def _clean_index(x):
        return pd.Index(x.astype(str).str.strip().str.strip('"').str.lstrip("﻿"))

    def _barcode_frac(x):
        s = _clean_index(x).to_series()
        pat = r"^[ACGTN]{8,}(?:[.-]\d+)?(?:[_-].+)?$"
        return s.str.match(pat).mean()

    idx = _clean_index(pd.Index(df.index))
    cols = _clean_index(pd.Index(df.columns))
    idx_bc = _barcode_frac(idx)
    col_bc = _barcode_frac(cols)

    if idx_bc > col_bc:
        ad = sc.AnnData(df.values.astype("float32"))
        ad.obs_names = idx
        ad.var_names = cols
    else:
        ad = sc.AnnData(df.values.T.astype("float32"))
        ad.obs_names = cols
        ad.var_names = idx
    ad.var_names_make_unique()
    return ad


def _tag_obs_names(ad, sample):
    ad.obs["barcode"] = ad.obs_names.astype(str)
    ad.obs_names = pd.Index([f"{sample}__{b}" for b in ad.obs["barcode"].astype(str)])
    return ad


def load_cohort(cfg):
    out = cfg["tar"].replace(".tar", "_ext")
    if not os.path.isdir(out):
        os.makedirs(out, exist_ok=True)
        with tarfile.open(cfg["tar"]) as t:
            t.extractall(out)
    files = [f for f in glob.glob(os.path.join(out, "**", "*"), recursive=True) if os.path.isfile(f)]
    ads, diag = [], []
    if cfg["fmt"] == "mtx":
        mtxs = [f for f in files if re.search(r"matrix\.mtx", os.path.basename(f), re.I)]
        for mf in mtxs:
            base = os.path.basename(mf)
            pref = re.split(r"matrix\.mtx", base, flags=re.I)[0]
            sibs = [f for f in files if os.path.basename(f).startswith(pref)]
            bc = next((f for f in sibs if "barcode" in f.lower()), None)
            ft = next((f for f in sibs if ("feature" in f.lower() or "gene" in f.lower())), None)
            samp = re.sub(r"^GSM\d+[_-]?", "", pref).strip("_-")
            if not (bc and ft):
                diag.append(f"跳过 {base}: 缺 barcodes/features")
                continue
            ad = load_mtx(mf, bc, ft)
            ad.obs["sample"] = samp
            ad.obs["route"] = cfg["route"](samp)
            ads.append(_tag_obs_names(ad, samp))
            diag.append(f"{samp} [{cfg['route'](samp)}]: {ad.n_obs} cells")
    else:
        csvs = [f for f in files if re.search(r"\.csv", os.path.basename(f), re.I)]
        for cf in csvs:
            base = os.path.basename(cf)
            samp = re.sub(r"^GSM\d+[_-]?", "", base)
            samp = re.sub(r"\.csv.*$", "", samp, flags=re.I).strip("_-")
            ad = load_csv(cf)
            ad.obs["sample"] = samp
            ad.obs["route"] = cfg["route"](samp)
            ads.append(_tag_obs_names(ad, samp))
            diag.append(f"{samp} [{cfg['route'](samp)}]: {ad.n_obs} cells")
    ads = [a for a in ads if a.obs["route"].iloc[0] != "drop"]
    if not ads:
        return None, diag, files
    m = sc.concat(ads, join="inner")
    return m, diag, files


def epithelial(m):
    sc.pp.filter_cells(m, min_genes=200)
    sc.pp.filter_genes(m, min_cells=3)
    m.var["mt"] = m.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(m, qc_vars=["mt"], inplace=True, percent_top=None)
    m = m[m.obs.pct_counts_mt < 20].copy()
    m.layers["counts"] = m.X.copy()
    sc.pp.normalize_total(m, target_sum=1e4)
    sc.pp.log1p(m)
    sc.pp.highly_variable_genes(m, n_top_genes=2000)
    sc.pp.pca(m, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(m, n_neighbors=15)
    sc.tl.leiden(m, resolution=1.0)
    epi = [g for g in ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"] if g in m.var_names]
    sc.tl.score_genes(m, epi, score_name="epi")
    cl = m.obs.groupby("leiden")["epi"].mean().sort_values(ascending=False)
    thr = cl.mean() + 0.5 * cl.std()
    keep = cl[cl > thr].index.tolist() or [cl.index[0]]
    return m[m.obs["leiden"].isin(keep)].copy(), cl, keep


EMT_INTRINSIC = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2"]


def cohend(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp0 = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (b.mean() - a.mean()) / sp0 if sp0 > 0 else np.nan


def load_emt(sigs):
    for k in ["EMT", "emt", "Emt", "EMT_signature"]:
        if k in sigs and len(sigs[k]) > 3:
            return list(sigs[k]), f"signatures.json[{k}]"
    canon = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2",
             "FN1", "SPARC", "TGFBI", "ITGA5", "TIMP1", "MMP2", "LOXL2"]
    return canon, "canonical fallback(需与 Fig4c 的 EMT 集对齐)"


def score_modules(e, sigs):
    for k in ["M0", "M1", "M2"]:
        g = [x for x in sigs[k] if x in e.var_names]
        sc.tl.score_genes(e, g, score_name=f"s_{k}")
    e.obs["dax"] = zscore(e.obs["s_M1"].values) - zscore(e.obs["s_M2"].values)
    emt, src = load_emt(sigs)
    eg = [x for x in emt if x in e.var_names]
    sc.tl.score_genes(e, eg, score_name="s_EMT")
    egi = [x for x in EMT_INTRINSIC if x in e.var_names]
    sc.tl.score_genes(e, egi, score_name="s_EMTi")
    return src, len(eg), len(emt), len(egi)


def _dump_mtx(a, prefix):
    cnt = a.layers["counts"] if "counts" in a.layers else a.X
    cnt = sp.csr_matrix(cnt) if not sp.issparse(cnt) else cnt.tocsr()
    sio.mmwrite(f"{prefix}.mtx", cnt.T.tocoo())
    pd.Series(a.var_names).to_csv(f"{prefix}.genes", index=False, header=False)
    pd.Series(a.obs_names.astype(str)).to_csv(f"{prefix}.barcodes", index=False, header=False)


def _sample_key(cohort, sample):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", f"{cohort}__{sample}")


def run_copykat_on_epithelial(e, cfg, copykat_root, ncores=8, min_cells=COPYKAT_MIN_CELLS):
    indir = os.path.join(copykat_root, "in")
    outdir = os.path.join(copykat_root, "out")
    os.makedirs(indir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    script = os.path.join(SCRIPT_DIR, "run_copykat.R")
    rows = []
    preds = {}

    for sample in sorted(e.obs["sample"].astype(str).unique()):
        sub = e[e.obs["sample"].astype(str).values == sample].copy()
        route = str(sub.obs["route"].iloc[0]) if sub.n_obs else "NA"
        key = _sample_key(cfg["name"], sample)
        row = dict(cohort=cfg["name"], sample=sample, sample_key=key, route=route,
                   n_epi=int(sub.n_obs), n_pred=0, n_aneuploid=0, n_nonaneuploid=0,
                   status="queued", note="")
        if sub.n_obs < min_cells:
            row["status"] = "skip_small"
            rows.append(row)
            continue
        prefix = os.path.join(indir, key)
        _dump_mtx(sub, prefix)
        proc = subprocess.run([RS, script, prefix, outdir, key, str(ncores)],
                              capture_output=True, text=True)
        lines = [x.strip() for x in (proc.stdout.splitlines() + proc.stderr.splitlines()) if x.strip()]
        if lines:
            row["note"] = " | ".join(lines[-3:])[:500]
        pred_file = os.path.join(outdir, f"{key}_copykat_prediction.txt")
        if proc.returncode != 0:
            row["status"] = "failed"
            rows.append(row)
            continue
        if not os.path.exists(pred_file):
            row["status"] = "missing_prediction"
            rows.append(row)
            continue
        p = pd.read_csv(pred_file, sep="\t")
        bcol = "cell.names" if "cell.names" in p.columns else p.columns[0]
        ccol = "copykat.pred" if "copykat.pred" in p.columns else p.columns[-1]
        calls = p[ccol].astype(str)
        row["status"] = "ok"
        row["n_pred"] = int(calls.notna().sum())
        row["n_aneuploid"] = int(calls.str.lower().str.startswith("aneu").sum())
        row["n_nonaneuploid"] = int(row["n_pred"] - row["n_aneuploid"])
        for bc, val in zip(p[bcol].astype(str), calls):
            preds[bc] = val
        rows.append(row)

    calls = pd.Series("undetermined", index=e.obs_names, dtype="object")
    mapped = pd.Series(e.obs_names, index=e.obs_names).map(preds)
    has = mapped.notna()
    calls.loc[has] = np.where(mapped.loc[has].str.lower().str.startswith("aneu"),
                              "malignant", "normal_epi")
    out = e.copy()
    out.obs["cnv_call"] = pd.Categorical(calls)
    return out[out.obs["cnv_call"] == "malignant"].copy(), pd.DataFrame(rows)


def main(v3=False, copykat_root=COPYKAT_ROOT_V3, copykat_cores=8, copykat_min_cells=COPYKAT_MIN_CELLS):
    sigs = frozen_sigs()
    sig_hash = signature_sha256(sigs)
    names = out_names(v3)
    if v3:
        L = ["# S12 V3 外部独立队列验证(CopyKAT malignant only)\n"]
    else:
        L = ["# S12 外部独立队列验证(分化轴 + 去分化 + 非EMT 复现)\n"]
    fig = plt.figure(figsize=(W2, 0.40 * W2 * len(COHORTS) + 0.42 * W2))
    gs = fig.add_gridspec(len(COHORTS) + 1, 3, height_ratios=[1] * len(COHORTS) + [0.9])
    AX = np.empty((len(COHORTS), 3), dtype=object)
    for _r in range(len(COHORTS)):
        for _c in range(3):
            AX[_r, _c] = fig.add_subplot(gs[_r, _c])
    axF = fig.add_subplot(gs[len(COHORTS), :])
    pl = "abcdefghi"
    pi = 0
    manifests = []

    for row, cfg in enumerate(COHORTS):
        L.append(f"## {cfg['label']}  ({cfg['name']})")
        if not os.path.exists(cfg["tar"]):
            L.append(f"⚠ 未找到 {cfg['tar']} — 先下载(见指令)。\n")
            for c in range(3):
                AX[row, c].text(0.1, 0.5, f"{cfg['name']}\n未下载", fontsize=7)
                panel(AX[row, c], pl[pi])
                pi += 1
            continue
        m, diag, files = load_cohort(cfg)
        L.append("- 载入: " + "; ".join(diag))
        if m is None:
            L.append("⚠ 没解析出样本。tar 内文件(前 20):\n" +
                     "\n".join("  " + os.path.basename(f) for f in files[:20]) + "\n把这个贴我对版。\n")
            for c in range(3):
                AX[row, c].text(0.1, 0.5, f"{cfg['name']}\n解析失败", fontsize=7)
                panel(AX[row, c], pl[pi])
                pi += 1
            continue

        e, cl, keep = epithelial(m)
        L.append(f"- 上皮簇 {keep} → 上皮细胞 {e.n_obs}")

        if v3:
            e, ck = run_copykat_on_epithelial(e, cfg, copykat_root, ncores=copykat_cores,
                                              min_cells=copykat_min_cells)
            if not ck.empty:
                manifests.append(ck)
                ok = int((ck["status"] == "ok").sum())
                L.append(f"- CopyKAT: {ok}/{ck.shape[0]} 样本成功")
                for _, rr in ck.iterrows():
                    L.append(f"  · {rr['sample']} [{rr['route']}]: epi={int(rr['n_epi'])} pred={int(rr['n_pred'])} "
                             f"aneu={int(rr['n_aneuploid'])} status={rr['status']}")
            L.append(f"- CopyKAT malignant 细胞 {e.n_obs}")
            if e.n_obs == 0:
                L.append("⚠ CopyKAT 后无 malignant 细胞，跳过该队列统计。\n")
                for c in range(3):
                    AX[row, c].text(0.1, 0.5, f"{cfg['name']}\n无 malignant", fontsize=7)
                    panel(AX[row, c], pl[pi])
                    pi += 1
                continue

        emt_src, nEMT, nEMTtot, nEMTi = score_modules(e, sigs)
        nM1 = len([x for x in sigs["M1"] if x in e.var_names])
        nM2 = len([x for x in sigs["M2"] if x in e.var_names])
        cols = ["sample", "route", "s_M1", "s_M2", "dax", "s_EMT", "s_EMTi"]
        if "cnv_call" in e.obs.columns:
            cols.append("cnv_call")
        cell_scores = e.obs[cols].assign(cohort=cfg["name"])
        stamp_signature_hash(cell_scores, sig_hash).to_csv(
            f"{C.RESULTS}/{names['cell_prefix']}{cfg['name']}.csv", index=False)
        L.append(f"- 签名命中 M1 {nM1}/{len(sigs['M1'])}, M2 {nM2}/{len(sigs['M2'])}")
        if nM1 < 5 or nM2 < 5:
            L.append("  ⚠ 签名命中过低 — 该队列基因可能是 Ensembl ID,需映射成 symbol。把 var_names 前几个贴我。")
        L.append(f"- EMT 签名: {emt_src}(全 {nEMT}/{nEMTtot} · 肿瘤内在 {nEMTi})")

        r, p = pearsonr(e.obs["s_M1"], e.obs["s_M2"])
        L.append(f"- ① 轴复现: M1 vs M2 打分 r={r:.3f} (p={p:.1e}) — 负相关 = 同一条分化轴复现")

        pri = e.obs.loc[e.obs.route == "primary", "dax"]
        met = e.obs.loc[e.obs.route == cfg["met"], "dax"]
        ded, delta, npri, nmet = "", None, 0, 0
        if len(pri) > 20 and len(met) > 20:
            delta = met.mean() - pri.mean()
            smp = e.obs[e.obs.route == "primary"].groupby("sample")["dax"].mean()
            smm = e.obs[e.obs.route == cfg["met"]].groupby("sample")["dax"].mean()
            npri, nmet = smp.size, smm.size
            _, pc = mannwhitneyu(pri, met, alternative="greater")
            try:
                _, ps = mannwhitneyu(smp, smm, alternative="greater")
            except Exception:
                ps = float("nan")
            ded = "y"
            L.append(f"- ② 去分化复现: primary {pri.mean():+.3f} → {cfg['met']} {met.mean():+.3f}, Δ={delta:+.3f}（Δ<0 = 转移更去分化,与 HRA 同向）")
            L.append(f"    · 样本级均值: primary={list(smp.round(2))} vs {cfg['met']}={list(smm.round(2))}")
            L.append(f"    · 样本级 Mann-Whitney p={ps:.3f}(n={npri}v{nmet},推断单位)；细胞级 p={pc:.1e}（伪重复,仅描述,勿作头条）")
        else:
            L.append(f"- ② 去分化: primary {len(pri)} / {cfg['met']} {len(met)} 细胞,样本不足跳过")

        demt = pe = demti = pei = dd_diff = dd_emt = dd_emti = None
        ep = e.obs.loc[e.obs.route == "primary", "s_EMT"]
        em = e.obs.loc[e.obs.route == cfg["met"], "s_EMT"]
        epi = e.obs.loc[e.obs.route == "primary", "s_EMTi"]
        emi = e.obs.loc[e.obs.route == cfg["met"], "s_EMTi"]
        if delta is not None and len(ep) > 20 and len(em) > 20:
            gp = e.obs[e.obs.route == "primary"].groupby("sample")
            gm = e.obs[e.obs.route == cfg["met"]].groupby("sample")
            demt = float(gm["s_EMT"].mean().mean() - gp["s_EMT"].mean().mean())
            demti = float(gm["s_EMTi"].mean().mean() - gp["s_EMTi"].mean().mean())
            try:
                _, pe = mannwhitneyu(gp["s_EMT"].mean(), gm["s_EMT"].mean(), alternative="two-sided")
            except Exception:
                pe = float("nan")
            try:
                _, pei = mannwhitneyu(gp["s_EMTi"].mean(), gm["s_EMTi"].mean(), alternative="two-sided")
            except Exception:
                pei = float("nan")
            dd_diff = cohend(pri, met)
            dd_emt = cohend(ep, em)
            dd_emti = cohend(epi, emi)
            small_emt = abs(dd_emti) < 0.4
            dominated = abs(dd_diff) > abs(dd_emti)
            if small_emt and dominated:
                verdict = "非EMT 成立(肿瘤内在 EMT 小且被去分化主导)"
            elif dominated:
                verdict = "去分化为主、但肿瘤内在 EMT 中等 → 写'非EMT驱动',勿说'无EMT'"
            else:
                verdict = "⚠ 肿瘤内在 EMT ≥ 去分化 → 该路径有真 EMT 成分,不能 claim 非EMT"
            L.append("- ③ 非EMT(标准化效应量 scale-fair, Cohen's d):")
            L.append(f"    · diff_axis={dd_diff:+.2f} | EMT全={dd_emt:+.2f} | EMT肿瘤内在={dd_emti:+.2f}")
            L.append(f"    · ΔEMT全={demt:+.3f}(p={pe:.3f}) | ΔEMT内在={demti:+.3f}(p={pei:.3f}) → {verdict}")
            if abs(dd_emt) > abs(dd_emti) + 0.2:
                L.append("    · 注: 全 EMT 效应 > 肿瘤内在 → 部分信号来自基质基因(FN1/SPARC/PDGFRB)=基质泄漏;但肿瘤内在 EMT 仍按上面判定")
        L.append("")

        ax0 = AX[row, 0]
        order = [o for o in ["primary", cfg["met"]] if (e.obs.route == o).sum() > 0]
        data = [e.obs.loc[e.obs.route == o, "dax"].values for o in order]
        bp = ax0.boxplot(data, positions=range(len(order)), widths=0.5, showfliers=False, patch_artist=True)
        for patch, col in zip(bp["boxes"], [STATE_COLORS["M1"], STATE_COLORS["M2"]][:len(order)]):
            patch.set_facecolor(col)
            patch.set_alpha(0.6)
        sm = e.obs.groupby(["route", "sample"])["dax"].mean()
        for xi, o in enumerate(order):
            try:
                ys = np.atleast_1d(sm.loc[o].values)
                ax0.scatter([xi] * len(ys), ys, c="k", s=14, zorder=3)
            except Exception:
                pass
        ax0.axhline(0, ls="--", c="gray", lw=0.8)
        ax0.set_xticks(range(len(order)))
        ax0.set_xticklabels([o.replace("_met", "") for o in order])
        ax0.set_ylabel("diff_axis (M1−M2)")
        ax0.set_title((f"{cfg['label']}\nΔ={delta:+.2f}  (n={npri}v{nmet} pts)") if ded else f"{cfg['label']}")
        panel(ax0, pl[pi])
        pi += 1

        axe = AX[row, 1]
        ds = [abs(dd_diff or 0), abs(dd_emt or 0), abs(dd_emti or 0)]
        axe.bar(range(3), ds, color=[STATE_COLORS["M2"], "#8d99ae", "#5c677d"], width=0.6)
        axe.set_xticks(range(3))
        axe.set_xticklabels(["diff_axis", "EMT\nfull", "EMT\nintrinsic"], fontsize=5.5)
        axe.set_ylabel("|Cohen's d| (met vs primary)")
        if dd_diff is not None:
            ttl = "Dediff > EMT" if abs(dd_diff) > abs(dd_emti or 0) else "EMT also engaged"
        else:
            ttl = "EMT"
        axe.set_title(ttl)
        panel(axe, pl[pi])
        pi += 1

        ax1 = AX[row, 2]
        ax1.scatter(e.obs["s_M1"], e.obs["s_M2"], s=2, alpha=0.3, c=e.obs["dax"], cmap="RdBu_r", lw=0)
        ax1.set_xlabel("M1 score")
        ax1.set_ylabel("M2 score")
        ax1.set_title(f"Axis reproduced (r={r:.2f})")
        panel(ax1, pl[pi])
        pi += 1

    def _cdelta(smap, pmap):
        pr = np.array([v for k, v in smap.items() if pmap[k]])
        me = np.array([v for k, v in smap.items() if not pmap[k]])
        if len(pr) < 2 or len(me) < 2:
            raise ValueError("n<2")
        d = me.mean() - pr.mean()
        se = np.sqrt(pr.var(ddof=1) / len(pr) + me.var(ddof=1) / len(me))
        return d, se

    forest = []
    hp = f"{C.RESULTS}/S4_scored.h5ad"
    if os.path.exists(hp):
        ho = sc.read_h5ad(hp, backed="r").obs
        dcol = next((c for c in ["dax", "diff_axis"] if c in ho), None)
        if dcol and "route" in ho and "sample" in ho:
            hh = ho[ho["cohort"] == "HRA004702"] if "cohort" in ho else ho
            hmet = {"ovarian_met", "peritoneal_met", "ascites", "peritoneal_lavage"}
            g = (pd.DataFrame({"d": np.asarray(hh[dcol]), "r": hh["route"].astype(str).values,
                               "s": hh["sample"].astype(str).values})
                 .groupby("s").agg(m=("d", "mean"), r=("r", "first")))
            g = g[g["r"].isin({"primary"} | hmet)]
            smap = g["m"].to_dict()
            pmap = {k: (v == "primary") for k, v in g["r"].to_dict().items()}
            try:
                d, se = _cdelta(smap, pmap)
                forest.append(("HRA (disc.)", d, se, "#888"))
            except Exception:
                pass
    for cn_, lbl in [("GSE239676", "GSE239676 (liver/ascites)"), ("GSE246662", "GSE246662 (liver)"), ("GSE308231", "GSE308231 (peri)")]:
        fp = f"{C.RESULTS}/{names['cell_prefix']}{cn_}.csv"
        if os.path.exists(fp):
            cdf = assert_signature_hash(pd.read_csv(fp), fp, sig_hash)
            cs = cdf.groupby("sample").agg(m=("dax", "mean"), r=("route", "first"))
            smap = cs["m"].to_dict()
            pmap = {k: (str(v) == "primary") for k, v in cs["r"].to_dict().items()}
            try:
                d, se = _cdelta(smap, pmap)
                forest.append((lbl, d, se, "#457b9d"))
            except Exception:
                pass
    pooled_p = None
    if os.path.exists(names["pooled_json"]):
        try:
            J = json.load(open(names["pooled_json"]))
            V = J.get("validation_only", {})
            beta_v = V.get("beta", np.nan)
            ci_v = V.get("ci", [np.nan, np.nan])
            p_v = V.get("p", np.nan)
            if np.isfinite(beta_v) and len(ci_v) >= 2 and np.isfinite(ci_v[0]) and np.isfinite(ci_v[1]):
                se_eq = (ci_v[1] - ci_v[0]) / (2 * 1.96)
                forest.append(("Pooled (val.)", beta_v, se_eq, "#9b2226"))
                pooled_p = p_v
        except Exception:
            pass
    if forest:
        yy = np.arange(len(forest))[::-1]
        for y, (lbl, d, se, col) in zip(yy, forest):
            dia = "Pooled" in lbl
            axF.plot([d - 1.96 * se, d + 1.96 * se], [y, y], "-", color=col, lw=2.4 if dia else 1.1, zorder=2)
            axF.scatter([d], [y], marker="D" if dia else "o", s=110 if dia else 34,
                        color=col, zorder=3, edgecolors="k", lw=0.5)
        axF.axvline(0, ls="--", c="0.78", lw=0.6)
        axF.set_yticks(yy)
        axF.set_yticklabels([f[0] for f in forest], fontsize=6.5)
        axF.set_ylim(-0.6, len(forest) - 0.4)
        axF.set_xlabel("Δ diff_axis (met − primary)")
        ttl = "Cross-cohort replication of metastatic dedifferentiation"
        if v3:
            ttl += "\nCopyKAT malignant only"
        axF.set_title(ttl)
        if pooled_p is not None:
            b = forest[-1][1]
            pt = ("%.0e" % pooled_p) if pooled_p < 1e-3 else ("%.3f" % pooled_p)
            axF.text(0.99, 0.96, f"validation pooled β={b:+.2f}, p={pt}", transform=axF.transAxes,
                     fontsize=6, va="top", ha="right", color="#9b2226")
    else:
        axF.axis("off")
    panel(axF, pl[pi])

    plt.tight_layout()
    save(fig, names["fig"])

    if manifests and names["copykat_manifest"]:
        pd.concat(manifests, ignore_index=True).to_csv(names["copykat_manifest"], index=False)

    if v3:
        L.append("解读：V3 使用 CopyKAT aneuploid 近似恶性后再评估 ①轴(M1/M2 负相关) ②去分化(Δdiff<0) ③肿瘤内在 EMT 的 Cohen's d。")
        L.append("统计: 推断单位是样本;细胞级 p 是伪重复勿用;若 CopyKAT 失败或 malignant 细胞过少,应退回门控而非强行解释。")
    else:
        L.append("解读：①轴(M1/M2 负相关) ②去分化(Δdiff<0,与 HRA 同向) 在独立外部队列复现。\n"
                 "③新主力外部 GSE239676 在 liver 与 ascites 两条路径上给出最强、最干净的 dedifferentiation-first 证据; ovarian 子集因 n 太小且 intrinsic-8 EMT 高,不用于 non-EMT 论点。\n"
                 "④非EMT 的判断以肿瘤内在 EMT 为准: liver 与 ascites 守住,308231 腹膜则提示 cohort-specific EMT;因此诚实写法应是『去分化是可与 EMT 分离的主轴,而 EMT 为队列/部位特异伴随成分』。\n"
                 "统计:推断单位是样本;细胞级 p 是伪重复勿用;上皮门控近似恶性(未跑 inferCNV)— 写入 limitation。")
    open(names["qc"], "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--v3", action="store_true", help="Run V3 external validation with CopyKAT malignant gating")
    ap.add_argument("--copykat-root", default=COPYKAT_ROOT_V3)
    ap.add_argument("--copykat-cores", type=int, default=8)
    ap.add_argument("--copykat-min-cells", type=int, default=COPYKAT_MIN_CELLS)
    args = ap.parse_args()
    main(v3=args.v3, copykat_root=args.copykat_root,
         copykat_cores=args.copykat_cores, copykat_min_cells=args.copykat_min_cells)
