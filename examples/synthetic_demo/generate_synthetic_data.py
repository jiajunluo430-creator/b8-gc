#!/usr/bin/env python
"""Generate a fully synthetic AnnData for the Stage-2 frozen-signature demo.

No real or patient data is used anywhere in this script. Cell counts are
drawn from a simple negative-binomial-like model, with each cell's "true"
latent programme (M0/M1/M2) upweighting its own frozen signature genes.
Route composition is skewed so that metastatic routes are enriched for M2
(low-differentiation) relative to primary tumour, purely to give the demo a
non-trivial, interpretable diff_axis pattern to recover.

Usage:
    python generate_synthetic_data.py [out.h5ad]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

SIGNATURES_JSON = Path(__file__).resolve().parent.parent.parent / "signatures" / "signatures.json"

COHORTS = ["SynCohortA", "SynCohortB", "SynCohortC"]
ROUTES = ["primary", "ovarian_met", "peritoneal_met", "ascites"]
# Route-dependent probability of each latent programme (rows sum to 1).
ROUTE_PROGRAMME_PROB = {
    "primary":         {"M0": 0.20, "M1": 0.55, "M2": 0.25},
    "ovarian_met":      {"M0": 0.15, "M1": 0.35, "M2": 0.50},
    "peritoneal_met":   {"M0": 0.15, "M1": 0.30, "M2": 0.55},
    "ascites":          {"M0": 0.10, "M1": 0.25, "M2": 0.65},
}
N_SAMPLES_PER_ROUTE = 4
N_CELLS_PER_SAMPLE = 60
N_FILLER_GENES = 120


def load_signature_genes():
    sigs = json.load(open(SIGNATURES_JSON))
    all_genes = sorted({g for genes in sigs.values() for g in genes})
    return sigs, all_genes


def build_gene_universe(sig_genes, rng):
    filler = [f"FILLER{i:04d}" for i in range(N_FILLER_GENES)]
    return sig_genes + filler


def simulate(rng, sigs, genes):
    gene_idx = {g: i for i, g in enumerate(genes)}
    rows = []
    obs = []
    sample_counter = 0
    for cohort in COHORTS:
        for route in ROUTES:
            for s in range(N_SAMPLES_PER_ROUTE):
                sample_counter += 1
                sample_id = f"{cohort}_{route}_S{s+1}"
                patient_id = f"{cohort}_P{s+1}"
                programme_probs = ROUTE_PROGRAMME_PROB[route]
                states = list(programme_probs.keys())
                probs = list(programme_probs.values())
                for _ in range(N_CELLS_PER_SAMPLE):
                    true_state = rng.choice(states, p=probs)
                    baseline = rng.negative_binomial(4, 0.35, size=len(genes)).astype(float)
                    boost_genes = rng.choice(sigs[true_state], size=12, replace=False)
                    for g in boost_genes:
                        baseline[gene_idx[g]] += rng.negative_binomial(6, 0.4)
                    rows.append(baseline)
                    obs.append({
                        "cohort": cohort,
                        "route": route,
                        "sample": sample_id,
                        "patient": patient_id,
                        "true_state": true_state,
                    })
    X = np.vstack(rows)
    obs_df = pd.DataFrame(obs)
    obs_df.index = [f"cell_{i:06d}" for i in range(len(obs_df))]
    var_df = pd.DataFrame(index=genes)
    return ad.AnnData(X=X, obs=obs_df, var=var_df)


def main():
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "synthetic_demo.h5ad"
    rng = np.random.default_rng(0)
    sigs, sig_genes = load_signature_genes()
    genes = build_gene_universe(sig_genes, rng)
    adata = simulate(rng, sigs, genes)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_path)
    print(f"Wrote synthetic demo AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes -> {out_path}")


if __name__ == "__main__":
    main()
