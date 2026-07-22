from pathlib import Path
import sys


PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE))

import config as C  # noqa: E402


def test_gse163558_suffix_routes():
    expected = {
        "GSM_sample_PT1": "primary",
        "GSM_sample_NT2": "adjacent_normal",
        "GSM_sample_LN1": "LN_met",
        "GSM_sample_O1": "ovarian_met",
        "GSM_sample_P1": "peritoneal_met",
        "GSM_sample_Li1": "liver_met",
    }
    assert {sample: C.route_for_sample("GSE163558", sample) for sample in expected} == expected


def test_gse270680_suffix_routes():
    expected = {
        "GSM_sample_T": "primary",
        "GSM_sample_N": "adjacent_normal",
        "GSM_sample_L": "LN_met",
        "GSM_sample_P": "blood",
    }
    assert {sample: C.route_for_sample("GSE270680", sample) for sample in expected} == expected


def test_gse183904_accession_routes_and_unknown_guard():
    assert C.route_for_sample("GSE183904", "GSM5573484.csv.gz") == "peritoneal_met"
    assert C.route_for_sample("GSE183904", "GSM5573502.csv.gz") == "adjacent_normal"
    assert C.route_for_sample("GSE183904", "GSM5573467.csv.gz") == "primary"
    assert C.route_for_sample("GSE183904", "unexpected_sample") is None


def test_gse228598_accession_routes_and_unknown_guard():
    assert C.route_for_sample("GSE228598", "GSM7133742_matrix.mtx.gz") == "peritoneal_lavage"
    assert C.route_for_sample("GSE228598", "GSM7133745_matrix.mtx.gz") == "ascites"
    assert C.route_for_sample("GSE228598", "unexpected_sample") is None


def test_mapped_cohorts_use_verified_route_assignment():
    for cohort in C.SAMPLE_MAPPED_COHORTS:
        assert C.COHORTS[cohort]["sample_route"] == "verified"
        assert C.COHORTS[cohort]["route"] is None
