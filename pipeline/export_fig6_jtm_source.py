#!/usr/bin/env python
"""Export canonical Figure 6 source tables from public pipeline outputs.

Inputs are produced by ``s12_external_validation.py``,
``s12_gse239676_validation.py``, and ``s13_mixed_model.py``. No internal
figure-source intermediates are required.
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

import config as C
from sig_utils import assert_signature_hash


OUT = Path(os.environ.get("B8GC_FIGURE_SOURCE_DATA", os.path.join(C.WORK, "figure_source_data")))
RESULTS = Path(C.RESULTS)

ROUTE_LABELS = {
    "primary": "primary",
    "liver_met": "liver",
    "ascites": "ascites",
    "peritoneal_met": "peritoneal",
    "ovarian_met": "ovarian",
}
SMALL_ROUTES = {
    "GSE308231": ("Peritoneal", "peritoneal_met"),
    "GSE246662": ("Liver", "liver_met"),
}
PROGRAMME_ROWS = {
    ("GSE308231", "peritoneal_met"): "GSE308231 | peritoneal",
    ("GSE246662", "liver_met"): "GSE246662 | liver",
    ("GSE239676", "liver_met"): "GSE239676 | liver",
    ("GSE239676", "ascites"): "GSE239676 | ascites",
}
FOREST_ROUTES = {
    "GSE308231": "Peritoneal",
    "GSE246662": "Liver",
    "GSE239676": "Liver/ascites",
}


def required(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"required pipeline output is missing: {path}")
    return path


def ci_from_samples(values):
    values = np.asarray(values, dtype=float).reshape(-1)
    if len(values) < 2:
        return np.nan, np.nan
    se = values.std(ddof=1) / np.sqrt(len(values))
    mean = values.mean()
    return float(mean - 1.96 * se), float(mean + 1.96 * se)


def sample_p(primary, metastasis):
    if len(primary) == 0 or len(metastasis) == 0:
        return np.nan
    return float(mannwhitneyu(primary, metastasis, alternative="two-sided").pvalue)


def paired_p(diffs):
    diffs = np.asarray(diffs, dtype=float)
    if len(diffs) < 2 or np.allclose(diffs, 0):
        return np.nan
    return float(wilcoxon(diffs).pvalue)


def load_cells(path, cohort):
    path = required(path)
    df = assert_signature_hash(pd.read_csv(path), str(path)).copy()
    required_cols = {"sample", "route", "s_M1", "s_M2", "dax", "s_EMT", "s_EMTi"}
    missing = sorted(required_cols - set(df.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns {missing}")
    df["cohort"] = cohort
    return df


def sample_table(cells):
    keys = ["cohort", "sample", "route"]
    agg = {"dax": "mean", "s_EMT": "mean", "s_EMTi": "mean"}
    if "patient" in cells.columns:
        agg["patient"] = "first"
    return cells.groupby(keys, observed=True, sort=False).agg(agg).reset_index()


def standard_cell_table(df):
    out = df.copy()
    out["cell_id"] = [f"{cohort}_{i:06d}" for i, cohort in enumerate(out["cohort"], start=1)]
    out["condition"] = np.where(out["route"].eq("primary"), "primary", "metastasis")
    return out[["cohort", "cell_id", "s_M1", "s_M2", "condition"]].rename(
        columns={"s_M1": "M1", "s_M2": "M2"}
    )


def delta_stat(samples, met_route):
    primary = samples.loc[samples["route"].eq("primary"), "dax"]
    metastasis = samples.loc[samples["route"].eq(met_route), "dax"]
    delta = float(metastasis.mean() - primary.mean())
    lo, hi = ci_from_samples(metastasis.to_numpy() - primary.mean())
    return delta, lo, hi, sample_p(primary, metastasis), len(primary), len(metastasis)


def export_axis(cells_by_cohort):
    axis = pd.concat([standard_cell_table(df) for df in cells_by_cohort.values()], ignore_index=True)
    axis.to_csv(OUT / "Fig6_B_axis_recovery.csv", index=False)
    correlations = [
        {"cohort": cohort, "r": float(df["s_M1"].corr(df["s_M2"]))}
        for cohort, df in cells_by_cohort.items()
    ]
    pd.DataFrame(correlations).to_csv(OUT / "Fig6_B_correlations.csv", index=False)


def export_small_cohorts(samples_by_cohort):
    frames, stats = [], []
    for cohort, (route_label, met_route) in SMALL_ROUTES.items():
        df = samples_by_cohort[cohort]
        subset = df[df["route"].isin(["primary", met_route])].copy()
        subset["group"] = np.where(subset["route"].eq("primary"), "primary", "metastasis")
        subset["route"] = route_label
        subset.rename(columns={"sample": "sample_id", "dax": "diff_axis"}, inplace=True)
        frames.append(subset[["cohort", "route", "group", "sample_id", "diff_axis"]])
        delta, lo, hi, p, n_primary, n_metastasis = delta_stat(df, met_route)
        stats.append({
            "cohort": cohort, "route": route_label, "estimate": delta, "lo": lo, "hi": hi,
            "p": p, "n_primary": n_primary, "n_metastasis": n_metastasis,
        })
    pd.concat(frames, ignore_index=True).to_csv(OUT / "Fig6_CD_small_cohort_samples.csv", index=False)
    pd.DataFrame(stats).to_csv(OUT / "Fig6_CD_small_cohort_stats.csv", index=False)


def export_route_validation(samples):
    route_samples = samples[samples["route"].isin(ROUTE_LABELS)].copy()
    route_samples["route"] = route_samples["route"].map(ROUTE_LABELS)
    route_samples.rename(columns={"sample": "sample_id", "dax": "diff_axis"}, inplace=True)
    route_samples[["route", "sample_id", "diff_axis"]].to_csv(
        OUT / "Fig6_E_239676_route_samples.csv", index=False
    )

    stats = []
    for route in ["liver", "ascites", "ovarian"]:
        primary = route_samples.loc[route_samples["route"].eq("primary"), "diff_axis"]
        metastasis = route_samples.loc[route_samples["route"].eq(route), "diff_axis"]
        estimate = float(metastasis.mean() - primary.mean()) if len(metastasis) else np.nan
        lo, hi = ci_from_samples(metastasis.to_numpy() - primary.mean()) if len(metastasis) >= 3 else (np.nan, np.nan)
        stats.append({
            "route": route, "estimate": estimate, "lo": lo, "hi": hi,
            "p": sample_p(primary, metastasis) if len(metastasis) >= 3 else np.nan,
            "note": "descriptive" if len(metastasis) < 3 else "",
        })
    pd.DataFrame(stats).to_csv(OUT / "Fig6_E_239676_route_stats.csv", index=False)


def export_paired(samples):
    if "patient" not in samples.columns:
        raise ValueError("GSE239676 sample table lacks patient identifiers required for paired analysis")
    long_rows, stats = [], []
    primary = samples[samples["route"].eq("primary")][["patient", "dax"]].rename(columns={"dax": "primary"})
    for raw_route, route in [("liver_met", "liver"), ("ascites", "ascites")]:
        metastasis = samples[samples["route"].eq(raw_route)][["patient", "dax"]].rename(columns={"dax": "metastasis"})
        paired = primary.merge(metastasis, on="patient", how="inner").sort_values("patient")
        for row in paired.itertuples(index=False):
            long_rows.extend([
                {"patient": row.patient, "route": route, "stage": "primary", "diff_axis": row.primary},
                {"patient": row.patient, "route": route, "stage": "metastasis", "diff_axis": row.metastasis},
            ])
        diffs = paired["metastasis"] - paired["primary"]
        lo, hi = ci_from_samples(diffs)
        stats.append({
            "route": route, "n_pairs": len(diffs), "estimate": float(diffs.mean()),
            "lo": lo, "hi": hi, "p": paired_p(diffs),
        })
    pd.DataFrame(long_rows).to_csv(OUT / "Fig6_F_239676_paired_samples.csv", index=False)
    pd.DataFrame(stats).to_csv(OUT / "Fig6_F_239676_paired_stats.csv", index=False)


def export_programmes(samples_by_cohort):
    rows = []
    for (cohort, met_route), cohort_route in PROGRAMME_ROWS.items():
        df = samples_by_cohort[cohort]
        primary = df[df["route"].eq("primary")]
        metastasis = df[df["route"].eq(met_route)]
        for programme, column in [
            ("Dedifferentiation", "dax"),
            ("Full EMT", "s_EMT"),
            ("Intrinsic EMT", "s_EMTi"),
        ]:
            effect = float(metastasis[column].mean() - primary[column].mean())
            lo, hi = ci_from_samples(metastasis[column].to_numpy() - primary[column].mean())
            rows.append({
                "cohort_route": cohort_route, "programme": programme,
                "effect_raw": effect, "lo_raw": lo, "hi_raw": hi,
            })
    pd.DataFrame(rows).to_csv(OUT / "Fig6_G_programme_effects.csv", index=False)


def export_forest_and_loco(pooled_path):
    with required(pooled_path).open(encoding="utf-8") as fh:
        pooled = json.load(fh)

    forest_rows = []
    by_cohort = {row["cohort"]: row for row in pooled["cohort_rows"]}
    for cohort in ["GSE308231", "GSE246662", "GSE239676"]:
        row = by_cohort[cohort]
        forest_rows.append({
            "cohort": cohort, "route": FOREST_ROUTES[cohort], "estimate": row["diff"],
            "lo": row["lo"], "hi": row["hi"],
            "n_label": f"{int(row['npri']) + int(row['nmet'])} samples",
            "p": np.nan, "pooled": False,
        })
    reference = pooled["validation_only"]
    n_validation = sum(int(by_cohort[c]["npri"]) + int(by_cohort[c]["nmet"]) for c in FOREST_ROUTES)
    forest_rows.append({
        "cohort": "validation-only pooled", "route": "Three held-out cohorts",
        "estimate": reference["beta"], "lo": reference["ci"][0], "hi": reference["ci"][1],
        "n_label": f"{n_validation} samples", "p": reference["p"], "pooled": True,
    })
    pd.DataFrame(forest_rows).to_csv(OUT / "Fig6_H_validation_forest.csv", index=False)

    loco_rows = [{
        "omitted": "None", "estimate": reference["beta"], "lo": reference["ci"][0],
        "hi": reference["ci"][1], "n_label": f"{n_validation} samples", "reference": True,
    }]
    for row in pooled["loco_rows"]:
        loco_rows.append({
            "omitted": row["omitted"], "estimate": row["beta"], "lo": row["ci"][0],
            "hi": row["ci"][1], "n_label": f"{int(row['n_samples'])} samples", "reference": False,
        })
    meta = pooled["validation_fixed_meta"]
    loco_rows.append({
        "omitted": "all validation", "estimate": meta["delta"], "lo": meta["ci"][0],
        "hi": meta["ci"][1], "n_label": f"{n_validation} samples", "reference": False,
    })
    pd.DataFrame(loco_rows).to_csv(OUT / "Fig6_I_loco.csv", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells_by_cohort = {
        "GSE308231": load_cells(RESULTS / "S12_cells_GSE308231.csv", "GSE308231"),
        "GSE246662": load_cells(RESULTS / "S12_cells_GSE246662.csv", "GSE246662"),
        "GSE239676": load_cells(RESULTS / "S12_239676_cells.csv", "GSE239676"),
    }
    samples_by_cohort = {cohort: sample_table(df) for cohort, df in cells_by_cohort.items()}
    export_axis(cells_by_cohort)
    export_small_cohorts(samples_by_cohort)
    export_route_validation(samples_by_cohort["GSE239676"])
    export_paired(samples_by_cohort["GSE239676"])
    export_programmes(samples_by_cohort)
    export_forest_and_loco(RESULTS / "S13_pooled.json")
    print("Wrote canonical Figure 6 source CSVs to", OUT)


if __name__ == "__main__":
    main()
