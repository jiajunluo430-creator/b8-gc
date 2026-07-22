#!/usr/bin/env python
"""s0_qdp.py — Stage 0: QC / Doublet / Preprocess（逐 cohort）。

Outputs under B8GC_WORK_ROOT: h5ad/<cohort>.h5ad, per-sample QC tables,
and qc/S0_summary.md.

用法:  python s0_qdp.py --cohort GSE163558
       python s0_qdp.py --all

关键纪律（见 CLAUDE.md）:
  - 全 cohort 用 config.QC_PARAMS 的**同一套阈值**；偏离须在 qc 摘要记录理由。
  - obs['route'] 必须写对（HRA004702 按样本后缀映射）。
  - 逐 cohort 处理、即时释放内存。
"""
import os, gc, argparse, glob, re, gzip
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import config as C

sc.settings.n_jobs = 16
sc.settings.verbosity = 1


# ----------------------------- loaders -----------------------------
def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _read_lines(path):
    with _open(path) as fh:
        return [ln.rstrip("\n") for ln in fh]


def _read_mtx(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return sio.mmread(fh).tocsr()
    return sio.mmread(path).tocsr()


def _features_to_symbols(feat_path):
    """features.tsv: 第2列=symbol(v3 三列) 或单列=symbol(v2)。取 symbol。"""
    out = []
    for ln in _read_lines(feat_path):
        parts = ln.split("\t")
        out.append(parts[1] if len(parts) >= 2 else parts[0])
    return out


def _build_triplet(mtx_f, bc_f, feat_f, sample):
    M = _read_mtx(mtx_f)                              # 10x 约定: genes × cells
    bcs = _read_lines(bc_f)
    genes = _features_to_symbols(feat_f)
    if M.shape == (len(genes), len(bcs)):
        X = M.T.tocsr()
    elif M.shape == (len(bcs), len(genes)):
        X = M.tocsr()
    else:
        raise RuntimeError(f"{sample}: mtx {M.shape} 与 genes({len(genes)})/barcodes({len(bcs)}) 不符")
    a = ad.AnnData(X.astype("float32"),
                   obs=pd.DataFrame(index=[str(b) for b in bcs]),
                   var=pd.DataFrame(index=[str(g) for g in genes]))
    a.var_names_make_unique()
    a.obs["sample"] = sample
    return a


def _load_mtx10x(path):
    """10x 矩阵，支持三种布局:
       (a) 每样本子目录含标准三件套(matrix.mtx.gz/features/barcodes)；
       (b) 扁平目录的 GSM 前缀三件套(GSMxxx_matrix.mtx.gz / _barcodes / _features|genes)；
       (c) 扁平目录单样本标准三件套。"""
    adatas = []
    # (a) 子目录标准三件套
    subs = [d for d in sorted(os.listdir(path)) if os.path.isdir(os.path.join(path, d))]
    for s in subs:
        try:
            a = sc.read_10x_mtx(os.path.join(path, s), var_names="gene_symbols", make_unique=True)
            a.obs["sample"] = s; adatas.append(a)
        except Exception:
            pass
    if adatas:
        return ad.concat(adatas, join="outer", index_unique="-")

    # (b/c) 扁平目录：按前缀分组三件套
    files = sorted(os.listdir(path))
    mtx_files = [f for f in files if re.search(r"matrix\.mtx(\.gz)?$", f)]
    if not mtx_files:
        raise RuntimeError(
            f"no 10x matrix under {path}（无子目录三件套，也无 *matrix.mtx[.gz]）。"
            f"请 ls 该目录确认布局后告知")
    for mf in mtx_files:
        prefix = re.sub(r"matrix\.mtx(\.gz)?$", "", mf).rstrip("_.")   # 'GSMxxx_' 或 ''

        def _pick(keys):
            for f in files:
                if (prefix == "" or f.startswith(prefix)) and any(k in f for k in keys) \
                        and "matrix" not in f:
                    return os.path.join(path, f)
            return None
        bc = _pick(["barcode"]); feat = _pick(["feature", "genes", "gene"])
        if not bc or not feat:
            raise RuntimeError(f"{mf}: 缺配套 barcodes/features（prefix='{prefix}'）→ ls 目录告知")
        sample = prefix or os.path.basename(path.rstrip("/"))
        adatas.append(_build_triplet(os.path.join(path, mf), bc, feat, sample))
    return ad.concat(adatas, join="outer", index_unique="-")


def _load_txt(path):
    """目录下逐样本矩阵(genes×cells)：.txt/.tsv(.gz)=tab，.csv(.gz)=comma。
       已核对 GSE134520 与 GSE183904_raw 均为 genes×cells → 读入后转置；用稀疏省内存。"""
    import scipy.sparse as sp
    files = sorted(glob.glob(f"{path}/*.txt.gz") + glob.glob(f"{path}/*.txt")
                   + glob.glob(f"{path}/*.tsv.gz") + glob.glob(f"{path}/*.tsv")
                   + glob.glob(f"{path}/*.csv.gz") + glob.glob(f"{path}/*.csv"))
    if not files:
        raise RuntimeError(f"no matrix files (.txt/.tsv/.csv[.gz]) under {path}")
    adatas = []
    for f in files:
        sep = "," if ".csv" in os.path.basename(f) else "\t"
        df = pd.read_csv(f, sep=sep, index_col=0)              # genes × cells
        a = ad.AnnData(sp.csr_matrix(df.T.values.astype("float32")),
                       obs=pd.DataFrame(index=df.columns.astype(str)),
                       var=pd.DataFrame(index=df.index.astype(str)))   # → cells × genes
        a.obs["sample"] = os.path.basename(f).split(".")[0]
        adatas.append(a)
        del df
    return ad.concat(adatas, join="outer", index_unique="-")


def _looks_like_counts(X):
    x = X[:50].toarray() if hasattr(X, "toarray") else np.asarray(X[:50])
    return np.allclose(x, np.round(x)) and float(x.min()) >= 0


def _load_from_h5ad(meta, cohort):
    """RData/RDS：先在 R 端转 h5ad（见 TASK_S0_QDP.md），本函数只读 h5ad。
       X 必须像 counts；否则找 counts 层；再否则：allow_noncounts=True 才放行(pass-through)。"""
    h5 = meta.get("h5ad_hint")
    assert h5 and os.path.exists(h5), (
        f"{cohort}({meta['fmt']}) 需先在 R 端转 h5ad 到 {h5}（见 TASK_S0_QDP.md）")
    a = sc.read_h5ad(h5)
    if not _looks_like_counts(a.X):
        if "counts" in a.layers:
            a.X = a.layers["counts"].copy()
        elif meta.get("allow_noncounts"):
            print(f"  [WARN] {cohort}: X 非 counts 且无 counts 层，但 allow_noncounts=True "
                  f"→ 按已归一化对象 pass-through（不重过滤/不重归一化）")
        else:
            raise ValueError(
                f"{cohort}: X 非 counts 且无 counts 层。先用增强版 export_rds.R 看诊断"
                f"(assays/layers/integer)；若 RNA counts 确不可恢复，再在 config 给该队列 allow_noncounts=True")
    if "sample" not in a.obs:
        for c in ("orig.ident", "sample_id", "Sample", "patient", "samples"):
            if c in a.obs:
                a.obs["sample"] = a.obs[c]; break
        else:
            a.obs["sample"] = cohort
    return a


def _load_cellranger_multi(path):
    """Load one Cell Ranger filtered-feature matrix directory per sample."""
    adatas = []
    for s in sorted(os.listdir(path)):
        mtx = os.path.join(path, s, "filtered_feature_bc_matrix")
        if not os.path.isdir(mtx):
            mtx = os.path.join(path, s)
        try:
            a = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=True)
        except Exception as e:
            print(f"  [skip sample] {s}: {e}"); continue
        a.obs["sample"] = s
        adatas.append(a)
    return ad.concat(adatas, join="outer", index_unique="-")


LOADERS = {"mtx10x": _load_mtx10x, "txt": _load_txt,
           "cellranger_multi": _load_cellranger_multi}
# fmt in ("rdata","rds") 走 _load_from_h5ad（先 R 端转 h5ad）


# ----------------------------- route -----------------------------
def _route_from_suffix(s):
    s = str(s).upper()
    for suf, r in C.HRA_SUFFIX2ROUTE.items():
        if s.endswith(suf) or f"_{suf}" in s:
            return r
    return "unknown"


def _normalize_route(v):
    """归一 route：支持 'OM'、'P01_OM'(取末段后缀)、'ovarian'等措辞 → C.ROUTES。"""
    s = str(v).upper().strip()
    if s in C.HRA_SUFFIX2ROUTE:
        return C.HRA_SUFFIX2ROUTE[s]
    tail = re.split(r"[_\-.]", s)[-1]                  # P01_OM -> OM
    if tail in C.HRA_SUFFIX2ROUTE:
        return C.HRA_SUFFIX2ROUTE[tail]
    keys = {"PRIMAR": "primary", "TUMOR": "primary",
            "OVAR": "ovarian_met", "PERITON": "peritoneal_met",
            "ASCIT": "ascites", "ADJAC": "adjacent_normal", "NORMAL": "adjacent_normal",
            "LYMPH": "LN_met"}
    for k, r in keys.items():
        if k in s:
            return r
    return "unknown"


def assign_route(adata, cohort, meta):
    sr = meta.get("sample_route")
    if sr == "verified":
        adata.obs["route"] = adata.obs["sample"].map(
            lambda sample: C.route_for_sample(cohort, sample) or "unknown"
        )
    elif sr == "metadata":
        col = meta.get("route_col")
        if col and col in adata.obs:
            adata.obs["route"] = adata.obs[col].map(_normalize_route)
        else:
            print(f"  [WARN] {cohort}: route_col '{col}' 不在 obs，回退样本名后缀")
            adata.obs["route"] = adata.obs["sample"].map(_route_from_suffix)
    elif sr == "suffix":
        adata.obs["route"] = adata.obs["sample"].map(_route_from_suffix)
    else:
        adata.obs["route"] = meta["route"]
    adata.obs["cohort"] = cohort
    n_unknown = int((adata.obs["route"] == "unknown").sum())
    if n_unknown:
        print(f"  [WARN] {cohort}: {n_unknown} cells route=unknown → 检查 route_col/样本命名")
    return adata


# ----------------------------- QC core -----------------------------
def qc_one_sample(a, p):
    """单样本 QC + doublet。返回过滤后 adata 与保留统计。"""
    n0 = a.n_obs
    a.var_names_make_unique()
    a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    sc.pp.filter_cells(a, min_genes=p["min_genes"])
    a = a[a.obs.n_genes_by_counts >= p["min_genes"]]
    a = a[a.obs.total_counts >= p["min_counts"]]
    a = a[a.obs.pct_counts_mt <= p["max_pct_mt"]].copy()
    # doublet（逐样本）
    n_dbl = 0
    if p["scrublet"] and a.n_obs > 50:
        try:
            sc.external.pp.scrublet(a, verbose=False)
            n_dbl = int(a.obs.predicted_doublet.sum())
            a = a[~a.obs.predicted_doublet].copy()
        except Exception as e:
            print(f"    scrublet skip: {e}")
    return a, dict(n_raw=n0, n_kept=a.n_obs, n_doublet=n_dbl,
                   pct_kept=round(100 * a.n_obs / max(n0, 1), 1))


def process_cohort(cohort):
    meta = C.COHORTS[cohort]; p = C.QC_PARAMS
    out = f"{C.H5AD}/{cohort}.h5ad"
    if os.path.exists(out):
        print(f"[{cohort}] exists → skip"); return
    print(f"[{cohort}] loading ({meta['fmt']}) ...")
    adata = (_load_from_h5ad(meta, cohort) if meta["fmt"] in ("rdata", "rds")
             else LOADERS[meta["fmt"]](meta["path"]))

    # 验证集/已预处理对象的特殊路径（meta['preprocessed']）
    #   "qc_done_needs_norm": 已 QC(含去污染/doublet)，但需按标准管线归一化（与 discovery 一致）
    #   "lognorm":            已是 log 归一化，原样保留
    prep = meta.get("preprocessed")
    if prep in ("qc_done_needs_norm", "lognorm"):
        adata = assign_route(adata, cohort, meta)
        if meta.get("group_map"):
            adata.obs["ln_status"] = adata.obs["sample"].map(meta["group_map"]).fillna("NA")
        if prep == "qc_done_needs_norm":
            adata.layers["counts"] = adata.X.copy()      # 去污染分数 counts（count 尺度）
            sc.pp.normalize_total(adata, target_sum=p["target_sum"]); sc.pp.log1p(adata)
            note = "已 QC(去污染/doublet)，按标准 normalize_total+log1p 归一化（与 discovery 一致）"
        else:
            adata.layers["lognorm_input"] = adata.X.copy()
            note = "已 log 归一化，pass-through"
        adata.raw = adata
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=p["n_hvg"])   # 默认 flavor，跑在 log 数据上
        except Exception:
            pass
        adata.write(out)
        with open(f"{C.QC}/S0_summary.md", "a") as fh:
            fh.write(f"\n## {cohort}（{note}）\n")
            fh.write(f"- cells: {adata.n_obs}; routes: {dict(adata.obs.route.value_counts())}\n")
            if "ln_status" in adata.obs:
                fh.write(f"- ln_status: {dict(adata.obs.ln_status.value_counts())}\n")
        print(f"[{cohort}] {prep} done: {adata.n_obs} cells")
        del adata; gc.collect()
        return

    sc.pp.filter_genes(adata, min_cells=p["min_cells"])

    # 逐样本 QC
    rows, kept = [], []
    for s, sub in [(s, adata[adata.obs["sample"] == s].copy())
                   for s in adata.obs["sample"].unique()]:
        sub, stat = qc_one_sample(sub, p); stat["sample"] = s
        rows.append(stat)
        if sub.n_obs > 0:
            kept.append(sub)
    adata = ad.concat(kept, join="outer", index_unique="-")
    adata = assign_route(adata, cohort, meta)
    if meta.get("group_map"):                      # 3v3: 注入 LN± 分组
        adata.obs["ln_status"] = adata.obs["sample"].map(meta["group_map"]).fillna("NA")
        miss = sorted(set(adata.obs["sample"].unique()) - set(meta["group_map"]))
        if miss:
            print(f"  [WARN] {cohort}: 样本不在 group_map（ln_status=NA）: {miss}")

    # 归一化 + HVG（batch-aware on raw counts）
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=p["target_sum"]); sc.pp.log1p(adata)
    adata.raw = adata
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=p["n_hvg"], flavor="seurat_v3",
                                    layer="counts", batch_key="sample")
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=p["n_hvg"])

    adata.write(out)
    qc_df = pd.DataFrame(rows)[["sample", "n_raw", "n_kept", "n_doublet", "pct_kept"]]
    qc_df.to_csv(f"{C.QC}/{cohort}_qc.csv", index=False)
    _append_summary(cohort, adata, qc_df)
    print(f"[{cohort}] done: {adata.n_obs} cells, routes={sorted(adata.obs.route.unique())}")
    del adata, kept; gc.collect()


def _append_summary(cohort, adata, qc_df):
    with open(f"{C.QC}/S0_summary.md", "a") as fh:
        fh.write(f"\n## {cohort}\n")
        fh.write(f"- cells kept: {adata.n_obs}（mean pct_kept {qc_df.pct_kept.mean():.1f}%）\n")
        fh.write(f"- routes: {dict(adata.obs.route.value_counts())}\n")
        lo = qc_df[qc_df.pct_kept < 50]
        if len(lo):
            fh.write(f"- ⚠️ 低保留率样本(<50%): {lo['sample'].tolist()} → 人工核查\n")


if __name__ == "__main__":
    import sys, traceback
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort"); ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    targets = list(C.COHORTS) if a.all else [a.cohort]
    failed = []
    for c in targets:
        try:
            process_cohort(c)
        except Exception as e:
            traceback.print_exc()
            print(f"[{c}] FAILED: {e}")
            failed.append(c)
    if failed:
        sys.exit(f"FAILED cohorts: {failed}")   # 非零退出码 → 外层 || break 生效
