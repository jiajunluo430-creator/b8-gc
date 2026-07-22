#!/usr/bin/env python
"""用法: python build_h5ad.py <indir> <out.h5ad>
读取 export_rds.R 产出的 counts.mtx / genes.txt / barcodes.txt / meta.csv，
组装成 **counts 在 X** 的 h5ad（s0_qdp 的 _load_from_h5ad 会据此自检并续跑）。
仅依赖 scanpy/anndata/scipy（b8gc 环境内）。"""
import sys
import anndata as ad
import pandas as pd
import scipy.io as sio

if len(sys.argv) != 3:
    sys.exit("usage: build_h5ad.py <indir> <out.h5ad>")
indir, out = sys.argv[1], sys.argv[2]

M = sio.mmread(f"{indir}/counts.mtx").tocsr()                 # genes × cells
genes = [l.strip() for l in open(f"{indir}/genes.txt")]
cells = [l.strip() for l in open(f"{indir}/barcodes.txt")]
meta = pd.read_csv(f"{indir}/meta.csv", index_col=0)
meta.index = meta.index.astype(str)

obs = meta.reindex([str(c) for c in cells])                  # 按 barcodes 顺序对齐
A = ad.AnnData(M.T.tocsr(), obs=obs,
               var=pd.DataFrame(index=[str(g) for g in genes]))  # cells × genes
assert A.n_obs == len(cells) and A.n_vars == len(genes)
A.write(out)
print(f"h5ad {A.shape} (cells×genes), X=counts -> {out}")
