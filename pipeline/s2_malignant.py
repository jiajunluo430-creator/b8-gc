#!/usr/bin/env python
"""s2_malignant.py — Stage 2：恶性上皮细胞鉴定（CopyKAT 逐样本，非跳过项）。

为什么非跳过：若直接拿"上皮 leiden 簇"当恶性，跨队列"可复现态"可能只是上皮污染/doublet/
正常上皮的假象。必须先把真·恶性(非整倍体)细胞挑出来，复现性才有意义。
precancer(GSE134520)的化生上皮也在这步和真恶性分开(多为二倍体)。

流程（两步，--step）：
  prep    : 从 S1_integrated 取上皮细胞 → 回各队列 h5ad 取 counts → 按样本导出矩阵给 CopyKAT
            （CopyKAT 必须逐样本跑，避免把样本间差异误当 CNV）
  collect : 汇总 CopyKAT 预测 → 标 malignant=aneuploid → 写 S2_malignant.h5ad

中间：run_copykat.R 逐样本跑（见同目录脚本 + 下方 shell 循环）。

barcode 对齐：S1 用 sc.concat(index_unique="-") 追加了 "-<batch>" 后缀；
  原始 barcode = s1_index.rsplit("-",1)[0]（只剥最后一段，10x 的 "-1" 不受影响）。
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import scanpy as sc
import config as C

S1 = f"{C.RESULTS}/S1_integrated.h5ad"
OUT = f"{C.RESULTS}/S2_malignant.h5ad"
CKDIR = os.environ.get("B8GC_COPYKAT_SCRATCH", f"{C.WORK}/tmp/copykat")  # 中间：每样本输入/输出
MIN_EPI = 100                          # CopyKAT 逐样本最少上皮数；不足则跳过(标 undetermined)
EPI_LABEL = "Epithelial"


def _orig_bc(idx):
    return pd.Index([s.rsplit("-", 1)[0] for s in idx])


def prep():
    os.makedirs(f"{CKDIR}/in", exist_ok=True)
    os.makedirs(f"{CKDIR}/out", exist_ok=True)
    s1 = sc.read_h5ad(S1)
    epi = s1.obs[s1.obs["cell_type"] == EPI_LABEL][["cohort", "sample"]].copy()
    epi["orig"] = _orig_bc(epi.index)
    print(f"S1 上皮细胞: {len(epi)}  跨 {epi['cohort'].nunique()} 队列")

    manifest = []
    for coh in sorted(epi["cohort"].unique()):
        sub = epi[epi["cohort"] == coh]
        h5 = f"{C.H5AD}/{coh}.h5ad"
        a = sc.read_h5ad(h5)
        if "counts" not in a.layers:
            print(f"  [WARN] {coh} 无 counts 层 → 跳过(CopyKAT 需 counts)"); continue
        # 用原始 barcode 对齐
        a = a[a.obs_names.isin(set(sub["orig"]))].copy()
        bc2sample = dict(zip(sub["orig"], sub["sample"]))
        a.obs["sample"] = a.obs_names.map(bc2sample)
        for s, asub in _by_sample(a):
            if asub.n_obs < MIN_EPI:
                manifest.append((coh, s, asub.n_obs, "skip_small")); continue
            _dump_mtx(asub, f"{CKDIR}/in/{s}")
            manifest.append((coh, s, asub.n_obs, "queued"))
        del a
    mdf = pd.DataFrame(manifest, columns=["cohort", "sample", "n_epi", "status"])
    mdf.to_csv(f"{CKDIR}/manifest.csv", index=False)
    print(mdf.to_string(index=False))
    print(f"\n→ 已导出 {(mdf.status=='queued').sum()} 个样本待 CopyKAT；"
          f"{(mdf.status=='skip_small').sum()} 个上皮过少跳过。")
    print(f"下一步：对 {CKDIR}/in/ 下每个样本跑 run_copykat.R（见 TASK 指令）。")


def _by_sample(a):
    for s in a.obs["sample"].dropna().unique():
        yield s, a[a.obs["sample"] == s].copy()


def _dump_mtx(a, prefix):
    cnt = a.layers["counts"]
    cnt = sp.csr_matrix(cnt) if not sp.issparse(cnt) else cnt.tocsr()
    sio.mmwrite(f"{prefix}.mtx", cnt.T.tocoo())        # 写成 genes × cells
    pd.Series(a.var_names).to_csv(f"{prefix}.genes", index=False, header=False)
    pd.Series(a.obs_names).to_csv(f"{prefix}.barcodes", index=False, header=False)


def collect():
    s1 = sc.read_h5ad(S1)
    mdf = pd.read_csv(f"{CKDIR}/manifest.csv")
    # 读所有 CopyKAT 预测：每文件两列 cell.names,copykat.pred
    preds = {}
    for s in mdf.loc[mdf.status == "queued", "sample"]:
        f = f"{CKDIR}/out/{s}_copykat_prediction.txt"
        if not os.path.exists(f):
            print(f"  [missing] {s} 预测缺失（CopyKAT 可能失败）"); continue
        p = pd.read_csv(f, sep="\t")
        col = "copykat.pred" if "copykat.pred" in p.columns else p.columns[-1]
        bc = "cell.names" if "cell.names" in p.columns else p.columns[0]
        for b, v in zip(p[bc], p[col]):
            preds[b] = v        # b = 原始 barcode

    epi_mask = s1.obs["cell_type"] == EPI_LABEL
    orig = _orig_bc(s1.obs_names)
    call = pd.Series("non_epithelial", index=s1.obs_names)
    call[epi_mask.values] = "undetermined"
    mapped = pd.Series(orig, index=s1.obs_names).map(preds)
    has = mapped.notna() & epi_mask.values
    call[has] = np.where(mapped[has].str.lower().str.startswith("aneu"), "malignant", "normal_epi")
    s1.obs["cnv_call"] = pd.Categorical(call)
    print(s1.obs["cnv_call"].value_counts())

    mal = s1[s1.obs["cnv_call"] == "malignant"].copy()
    mal.write(OUT)
    with open(f"{C.QC}/S2_summary.md", "w") as fh:
        fh.write("# S2 Malignant Identification (CopyKAT)\n\n")
        fh.write(f"- 上皮总数: {int(epi_mask.sum())}; 恶性(aneuploid): {mal.n_obs}; "
                 f"正常上皮: {int((call=='normal_epi').sum())}; "
                 f"未定(样本过小/失败): {int((call=='undetermined').sum())}\n\n")
        fh.write("## malignant × cohort\n\n")
        fh.write(pd.crosstab(mal.obs.cohort, [1]*mal.n_obs).to_markdown() + "\n\n")
        fh.write("## malignant × route（这才是 Fig4 真正可用的跨路径恶性覆盖）\n\n")
        fh.write(mal.obs["route"].value_counts().to_markdown() + "\n")
    print("S2 完成 →", OUT)
    print("malignant × route:\n", mal.obs["route"].value_counts())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", choices=["prep", "collect"], required=True)
    a = ap.parse_args()
    prep() if a.step == "prep" else collect()
