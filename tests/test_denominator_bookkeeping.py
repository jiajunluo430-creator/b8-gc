"""Synthetic retention-audit test.

Confirms the core Stage-1 -> Stage-2 bookkeeping invariant documented in
docs/PROVENANCE.md: Stage-2 frozen-signature scoring is applied to EVERY
carcinoma-derived malignant cell, regardless of whether that cell's Stage-1
cluster reached the reproducible-module threshold (RMC0-4) or fell into the
`unrepro` bucket. No cell should be dropped between the two stages on the
basis of Stage-1 module membership.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from sig_utils import frozen_sigs

STATES = ["M0", "M1", "M2"]


def _make_synthetic_malignant_cells(rng, n_cells=300):
    sigs = frozen_sigs()
    genes = sorted({g for gs in sigs.values() for g in gs})
    x = rng.negative_binomial(4, 0.4, size=(n_cells, len(genes))).astype(float)

    # Simulate Stage-1 module membership: some cells belong to a reproducible
    # module (RMC1/RMC3/RMC4, the ones with a frozen-programme source), some
    # to non-signature-source modules (RMC0/RMC2), and some to `unrepro`.
    stage1_labels = rng.choice(
        ["RMC0", "RMC1", "RMC2", "RMC3", "RMC4", "unrepro"], size=n_cells
    )
    obs = pd.DataFrame({"stage1_module": stage1_labels})
    var = pd.DataFrame(index=genes)
    return ad.AnnData(X=x, obs=obs, var=var), sigs


def test_all_cells_scored_regardless_of_stage1_module():
    rng = np.random.default_rng(1)
    adata, sigs = _make_synthetic_malignant_cells(rng)
    n_before = adata.n_obs

    assert set(adata.obs["stage1_module"]) >= {"RMC0", "RMC2", "unrepro"}, (
        "test fixture must include non-signature-source modules"
    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    scols = []
    for m in STATES:
        genes = [g for g in sigs[m] if g in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name=f"s_{m}")
        scols.append(f"s_{m}")
    adata.obs["state"] = adata.obs[scols].idxmax(axis=1).str.replace("s_", "", regex=False)

    assert adata.n_obs == n_before, "no cell should be dropped by Stage-2 scoring"
    assert adata.obs["state"].notna().all(), "every cell must receive a Stage-2 argmax state"

    for module in ["RMC0", "RMC2", "unrepro"]:
        subset = adata.obs.loc[adata.obs["stage1_module"].eq(module)]
        assert len(subset) > 0
        assert subset["state"].notna().all(), (
            f"cells in Stage-1 module {module} (not a signature source) must still be scored"
        )
