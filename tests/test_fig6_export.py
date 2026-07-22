import json
from pathlib import Path

import pandas as pd

import export_fig6_jtm_source as ex


def _cells(cohort, met_route, with_patient=False):
    rows = []
    routes = ["primary", "primary", met_route, met_route]
    for i, route in enumerate(routes):
        row = {
            "cohort": cohort,
            "sample": f"P{i // 2}_{'P' if route == 'primary' else 'M'}_{i}",
            "route": route,
            "s_M1": 1.0 - i * 0.2,
            "s_M2": i * 0.2,
            "dax": 1.0 - i * 0.5,
            "s_EMT": i * 0.1,
            "s_EMTi": i * 0.05,
        }
        if with_patient:
            row["patient"] = f"P{i % 2}"
        rows.append(row)
    return pd.DataFrame(rows)


def test_exporter_uses_live_stage12_and_stage13_schemas(tmp_path, monkeypatch):
    monkeypatch.setattr(ex, "OUT", tmp_path)
    cells = {
        "GSE308231": _cells("GSE308231", "peritoneal_met"),
        "GSE246662": _cells("GSE246662", "liver_met"),
        "GSE239676": pd.concat(
            [
                _cells("GSE239676", "liver_met", with_patient=True),
                _cells("GSE239676", "ascites", with_patient=True).assign(
                    sample=lambda d: d["sample"] + "_A"
                ),
            ],
            ignore_index=True,
        ),
    }
    samples = {cohort: ex.sample_table(df) for cohort, df in cells.items()}
    ex.export_axis(cells)
    ex.export_small_cohorts(samples)
    ex.export_route_validation(samples["GSE239676"])
    ex.export_paired(samples["GSE239676"])
    ex.export_programmes(samples)

    pooled = {
        "validation_only": {"beta": -1.0, "ci": [-1.5, -0.5], "p": 0.001},
        "validation_fixed_meta": {"delta": -1.1, "ci": [-1.4, -0.8], "p": 0.0001},
        "cohort_rows": [
            {"cohort": c, "diff": -1.0, "lo": -2.0, "hi": 0.0, "npri": 2, "nmet": 2}
            for c in ["GSE308231", "GSE246662", "GSE239676"]
        ],
        "loco_rows": [
            {"omitted": c, "beta": -1.0, "ci": [-1.6, -0.4], "p": 0.01,
             "n_samples": 8, "retained_cohorts": []}
            for c in ["GSE308231", "GSE246662", "GSE239676"]
        ],
    }
    pooled_path = tmp_path / "S13_pooled.json"
    pooled_path.write_text(json.dumps(pooled), encoding="utf-8")
    ex.export_forest_and_loco(pooled_path)

    expected = {
        "Fig6_B_axis_recovery.csv", "Fig6_B_correlations.csv",
        "Fig6_CD_small_cohort_samples.csv", "Fig6_CD_small_cohort_stats.csv",
        "Fig6_E_239676_route_samples.csv", "Fig6_E_239676_route_stats.csv",
        "Fig6_F_239676_paired_samples.csv", "Fig6_F_239676_paired_stats.csv",
        "Fig6_G_programme_effects.csv", "Fig6_H_validation_forest.csv", "Fig6_I_loco.csv",
    }
    assert expected <= {path.name for path in tmp_path.glob("Fig6_*.csv")}
    assert list(pd.read_csv(tmp_path / "Fig6_I_loco.csv").columns) == [
        "omitted", "estimate", "lo", "hi", "n_label", "reference"
    ]


def test_exporter_has_no_legacy_figure_input_dependency():
    source = Path(ex.__file__).read_text(encoding="utf-8")
    for retired in [
        "Fig6_A_308231_diffaxis.csv", "Fig6_D_246662_diffaxis.csv",
        "Fig6_G_239676_route.csv", "Fig6_I_239676_paired.csv",
        "Fig6_J_forest.csv", "Fig6_L_loco.csv",
    ]:
        assert retired not in source
