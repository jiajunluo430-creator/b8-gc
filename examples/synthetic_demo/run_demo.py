#!/usr/bin/env python
"""End-to-end demo of Stage-2 frozen-signature scoring on synthetic data.

Mirrors the logic in pipeline/s4_crossroute.py: normalize -> per-programme
sc.tl.score_genes -> argmax programme assignment -> diff_axis = s_M1 - s_M2
-> route-level summary. Uses only the synthetic AnnData produced by
generate_synthetic_data.py; no real cohort data is read.

Usage:
    python generate_synthetic_data.py
    python run_demo.py
"""
import sys
from pathlib import Path

import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "pipeline"))
from sig_utils import frozen_sigs  # noqa: E402

STATES = ["M0", "M1", "M2"]


def main():
    h5ad_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "synthetic_demo.h5ad"
    if not h5ad_path.exists():
        raise SystemExit(f"{h5ad_path} not found — run generate_synthetic_data.py first.")

    adata = sc.read_h5ad(h5ad_path)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sigs = frozen_sigs()
    scols = []
    for m in STATES:
        genes = [g for g in sigs[m] if g in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name=f"s_{m}")
        scols.append(f"s_{m}")

    adata.obs["state"] = adata.obs[scols].idxmax(axis=1).str.replace("s_", "", regex=False)
    adata.obs["diff_axis"] = adata.obs["s_M1"] - adata.obs["s_M2"]

    print("=== Argmax programme composition by route ===")
    comp = adata.obs.groupby("route", observed=True)["state"].value_counts(normalize=True).unstack().round(3)
    print(comp.reindex(columns=STATES).to_string())

    print("\n=== Mean diff_axis by route (expect: primary highest, ascites lowest) ===")
    axis_by_route = adata.obs.groupby("route", observed=True)["diff_axis"].agg(["mean", "median", "count"]).round(3)
    print(axis_by_route.to_string())

    print("\n=== Recovery check: argmax state vs. simulated ground truth ===")
    agreement = (adata.obs["state"] == adata.obs["true_state"]).mean()
    print(f"argmax(state) == true_state in {agreement:.1%} of cells "
          "(not expected to be 100% — signature scoring is noisy by design; "
          "this just demonstrates the scoring pipeline runs end to end).")

    out_csv = h5ad_path.parent / "synthetic_demo_summary.csv"
    adata.obs[["cohort", "route", "sample", "patient", "true_state", "state", "diff_axis"]].to_csv(out_csv)
    print(f"\nWrote per-cell summary table -> {out_csv}")


if __name__ == "__main__":
    main()
