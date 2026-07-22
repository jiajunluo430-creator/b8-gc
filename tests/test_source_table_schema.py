"""Schema smoke tests for the signature table and scored-cell source tables."""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from sig_utils import (
    SIGNATURE_HASH_COL,
    SIGNATURE_STATES,
    SIGNATURE_TABLE,
    assert_signature_hash,
    frozen_sigs,
    load_supp_table_s2,
    signature_sha256,
    stamp_signature_hash,
)

STATES = ["M0", "M1", "M2"]


def test_supp_table_s2_schema():
    df = pd.read_csv(SIGNATURE_TABLE)
    assert {"state", "rank", "gene"} <= set(df.columns)
    assert set(df["state"].unique()) == set(SIGNATURE_STATES)


def test_load_supp_table_s2_round_trips_to_frozen_sigs():
    table_sigs = load_supp_table_s2()
    json_sigs = frozen_sigs()
    assert table_sigs == json_sigs


def test_scored_cell_table_schema():
    rng = np.random.default_rng(2)
    sigs = frozen_sigs()
    genes = sorted({g for gs in sigs.values() for g in gs})
    x = rng.negative_binomial(4, 0.4, size=(50, len(genes))).astype(float)
    adata = ad.AnnData(X=x, obs=pd.DataFrame(index=[f"c{i}" for i in range(50)]),
                        var=pd.DataFrame(index=genes))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    scols = []
    for m in STATES:
        gs = [g for g in sigs[m] if g in adata.var_names]
        sc.tl.score_genes(adata, gs, score_name=f"s_{m}")
        scols.append(f"s_{m}")
    adata.obs["state"] = adata.obs[scols].idxmax(axis=1).str.replace("s_", "", regex=False)
    adata.obs["diff_axis"] = adata.obs["s_M1"] - adata.obs["s_M2"]

    expected_cols = {"s_M0", "s_M1", "s_M2", "state", "diff_axis"}
    assert expected_cols <= set(adata.obs.columns)
    assert set(adata.obs["state"].unique()) <= set(STATES)

    sha = signature_sha256(sigs)
    stamped = stamp_signature_hash(adata.obs, sha)
    assert SIGNATURE_HASH_COL in stamped.columns
    assert_signature_hash(stamped, "synthetic-test-frame", expected_hash=sha)
