"""Shared frozen signature helpers.

[release note] SIGNATURE_TABLE and the signatures.json path read by
frozen_sigs() resolve, by default, to the signatures/ directory shipped
alongside pipeline/ in this release (overridable via B8GC_SIGNATURES_DIR),
rather than the original private figures/ export path.
"""
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

import config as C

EMT_INTRINSIC8 = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2"]
EMT_BROAD15 = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2",
               "FN1", "SPARC", "TGFBI", "ITGA5", "TIMP1", "MMP2", "LOXL2"]
EXPECTED_SIGNATURE_SHA256 = "b522d543b269a86476499ee2955c35a90cdcd3bce68c9d3f5f9c3634f7377324"
SIGNATURE_HASH_COL = "signature_sha256"
SIGNATURE_STATES = ["M0", "M1", "M2"]
SIGNATURES_DIR = Path(os.environ.get(
    "B8GC_SIGNATURES_DIR", Path(__file__).resolve().parent.parent / "signatures"))
SIGNATURE_TABLE = SIGNATURES_DIR / "SuppTable_S2_signatures.csv"


def canonical_signature_payload(sigs):
    return json.dumps(sigs, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def signature_sha256(sigs):
    return hashlib.sha256(canonical_signature_payload(sigs).encode()).hexdigest()


def validate_frozen_sigs(sigs, source="signatures.json"):
    sha = signature_sha256(sigs)
    if sha != EXPECTED_SIGNATURE_SHA256:
        raise RuntimeError(
            f"Frozen signature hash mismatch for {source}: observed {sha}, "
            f"expected {EXPECTED_SIGNATURE_SHA256}. Regenerate from canonical SuppTable S2 "
            "or intentionally update EXPECTED_SIGNATURE_SHA256."
        )
    return sha


def load_supp_table_s2(path=SIGNATURE_TABLE):
    df = pd.read_csv(path)
    required = {"state", "rank", "gene"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    df = df.loc[:, ["state", "rank", "gene"]].copy()
    df["state"] = df["state"].astype(str)
    df["rank"] = df["rank"].astype(int)
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()

    observed_states = list(df["state"].drop_duplicates())
    if observed_states != SIGNATURE_STATES:
        raise ValueError(f"Expected states {SIGNATURE_STATES}, observed {observed_states}")

    sigs = {}
    for state, sub in df.groupby("state", sort=False):
        sub = sub.sort_values("rank")
        ranks = sub["rank"].tolist()
        genes = sub["gene"].tolist()
        if ranks != list(range(1, 31)):
            raise ValueError(f"{state} ranks must be exactly 1..30, observed {ranks}")
        if len(set(genes)) != len(genes):
            dup = sorted({g for g in genes if genes.count(g) > 1})
            raise ValueError(f"{state} contains duplicated genes: {dup}")
        sigs[state] = genes
    return sigs


def assert_sigs_match_table_s2(sigs, source="signatures.json", table=SIGNATURE_TABLE):
    table_sigs = load_supp_table_s2(table)
    if sigs != table_sigs:
        raise RuntimeError(
            f"Frozen signatures in {source} do not match canonical Supplementary Table S2 at {table}. "
            "Run b8_pipeline/freeze_signatures.py and rerun downstream scorers."
        )
    return sigs


def frozen_sigs():
    path = os.environ.get("B8GC_SIGNATURES_JSON", str(SIGNATURES_DIR / "signatures.json"))
    sigs = json.load(open(path))
    validate_frozen_sigs(sigs, path)
    assert_sigs_match_table_s2(sigs, path)
    return sigs


def stamp_signature_hash(frame, sig_hash):
    out = frame.copy()
    out[SIGNATURE_HASH_COL] = sig_hash
    return out


def assert_signature_hash(frame, path, expected_hash=EXPECTED_SIGNATURE_SHA256):
    if SIGNATURE_HASH_COL not in frame.columns:
        raise RuntimeError(f"{path} lacks {SIGNATURE_HASH_COL}; rerun the S12 scorer with the current frozen signatures.")
    observed = set(frame[SIGNATURE_HASH_COL].dropna().astype(str))
    if observed != {expected_hash}:
        raise RuntimeError(
            f"{path} has stale signature hash {sorted(observed)}; expected {expected_hash}. "
            "Rerun the S12 scorer before pooling."
        )
    return frame
