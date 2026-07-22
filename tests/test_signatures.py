"""Smoke tests for the frozen M0/M1/M2 signature lock.

Regenerating or hand-editing signatures/signatures.json without updating
EXPECTED_SIGNATURE_SHA256 in pipeline/sig_utils.py should fail these tests.
"""
from sig_utils import (
    EXPECTED_SIGNATURE_SHA256,
    SIGNATURE_STATES,
    frozen_sigs,
    signature_sha256,
)


def test_frozen_sigs_loads_three_states():
    sigs = frozen_sigs()
    assert list(sigs.keys()) == SIGNATURE_STATES


def test_each_state_has_exactly_30_genes():
    sigs = frozen_sigs()
    for state, genes in sigs.items():
        assert len(genes) == 30, f"{state} has {len(genes)} genes, expected 30"


def test_no_duplicate_genes_within_a_state():
    sigs = frozen_sigs()
    for state, genes in sigs.items():
        assert len(set(genes)) == len(genes), f"{state} contains duplicate genes"


def test_canonical_hash_matches_lock():
    sigs = frozen_sigs()
    assert signature_sha256(sigs) == EXPECTED_SIGNATURE_SHA256


def test_frozen_sigs_raises_on_tampered_payload():
    import sig_utils

    sigs = frozen_sigs()
    tampered = {**sigs, "M0": ["NOTAREALGENE"] + sigs["M0"][1:]}
    try:
        sig_utils.validate_frozen_sigs(tampered, source="tampered-test")
    except RuntimeError:
        pass
    else:
        raise AssertionError("validate_frozen_sigs should reject a tampered signature payload")
