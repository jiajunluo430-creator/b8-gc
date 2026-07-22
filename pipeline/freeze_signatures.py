#!/usr/bin/env python
"""Freeze canonical malignant-state signatures from Supplementary Table S2."""
import json
from pathlib import Path

import config as C
from sig_utils import SIGNATURE_STATES as STATES, SIGNATURE_TABLE as TABLE, load_supp_table_s2, signature_sha256, validate_frozen_sigs

OUT_JSON = Path(C.RESULTS) / "signatures.json"
OUT_QC = Path(C.QC) / "SIG_freeze.md"


def main():
    sigs = load_supp_table_s2()
    sha = signature_sha256(sigs)
    validate_frozen_sigs(sigs, str(TABLE))
    OUT_JSON.write_text(json.dumps(sigs, indent=1) + "\n")

    lines = [
        "# 签名冻结（canonical Supplementary Table S2）",
        "",
        f"- Source: `{TABLE}`",
        f"- Output: `{OUT_JSON}`",
        f"- Canonical SHA256: `{sha}`",
        "",
    ]
    for state in STATES:
        lines += [f"## {state}", ", ".join(sigs[state]), ""]
    OUT_QC.write_text("\n".join(lines).rstrip() + "\n")

    print(f"Frozen signatures from {TABLE}")
    print(f"Wrote {OUT_JSON}")
    print(f"Canonical SHA256: {sha}")


if __name__ == "__main__":
    main()
