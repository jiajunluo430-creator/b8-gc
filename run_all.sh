#!/usr/bin/env bash
# Canonical re-run entrypoint for the stages included in this public release.
#
# This is a TRIMMED subset of the full internal pipeline (see docs/PROVENANCE.md
# and README.md for the excluded stages and why). It covers: QC + atlas
# construction, CopyKAT malignant identification, independent per-cohort
# clustering + MetaNeighbor raw-module construction, frozen-signature
# projection + differentiation-axis scoring, route/patient-paired summaries,
# EMT scoring, spatial Moran analysis, held-out validation + pooled/LOCO
# models, bulk negative-portability survival analyses, and Figure 6
# source-table export.
#
# All paths are resolved from B8GC_DATA_ROOT (raw input cohorts) and
# B8GC_WORK_ROOT (h5ad/results/qc/figures outputs) — set these before running,
# or rely on the ./data and ./work defaults. See README.md "Expected input
# formats" for what must exist under B8GC_DATA_ROOT before running S0.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="$ROOT/pipeline"
export PYTHONPATH="$PIPELINE${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg

export B8GC_DATA_ROOT="${B8GC_DATA_ROOT:-$ROOT/data}"
export B8GC_WORK_ROOT="${B8GC_WORK_ROOT:-$ROOT/work}"
export B8GC_SIGNATURES_DIR="${B8GC_SIGNATURES_DIR:-$ROOT/signatures}"

PYTHON_BIN="${PYTHON_BIN:-python}"
RSCRIPT_BIN="${B8GC_RSCRIPT:-${RSCRIPT_BIN:-Rscript}}"

run_py() {
  "$PYTHON_BIN" "$PIPELINE/$1" "${@:2}"
}

run_r() {
  "$RSCRIPT_BIN" "$PIPELINE/$1" "${@:2}"
}

# ---- Stage 0: raw RDS -> h5ad conversion for the two RDS-format cohorts ----
run_r export_rds.R \
  "$B8GC_DATA_ROOT/HRA004702/GCscrRNAlist.rds" \
  "$B8GC_WORK_ROOT/tmp/hra_export"
run_py build_h5ad.py "$B8GC_WORK_ROOT/tmp/hra_export" "$B8GC_WORK_ROOT/h5ad/HRA004702.raw.h5ad"

INCLUDE_INHOUSE=0
if [[ -f "$B8GC_DATA_ROOT/GC/GC_qc.rds" ]]; then
  run_r export_rds.R \
    "$B8GC_DATA_ROOT/GC/GC_qc.rds" \
    "$B8GC_WORK_ROOT/tmp/3v3_export"
  run_py build_h5ad.py "$B8GC_WORK_ROOT/tmp/3v3_export" "$B8GC_WORK_ROOT/h5ad/INHOUSE_3v3.raw.h5ad"
  INCLUDE_INHOUSE=1
else
  echo "[optional] INHOUSE_3v3 input not present; skipping its conversion and analysis"
fi

# ---- Stage 0: unified QC, per cohort ----
cohorts=(GSE163558 GSE270680 GSE228598 GSE134520 GSE183904 HRA004702)
if [[ "$INCLUDE_INHOUSE" == 1 ]]; then
  cohorts+=(INHOUSE_3v3)
fi
for cohort in "${cohorts[@]}"; do
  run_py s0_qdp.py --cohort "$cohort"
done

# ---- Stage 1: atlas integration ----
run_py s1_integrate.py

# ---- Stage 2: per-sample CopyKAT malignant identification ----
run_py s2_malignant.py --step prep
CKDIR="${B8GC_COPYKAT_SCRATCH:-$B8GC_WORK_ROOT/tmp/copykat}"
while IFS= read -r sample; do
  pred="$CKDIR/out/${sample}_copykat_prediction.txt"
  if [[ -s "$pred" ]]; then
    echo "[CopyKAT] $sample exists -> skip"
    continue
  fi
  run_r run_copykat.R "$CKDIR/in/$sample" "$CKDIR/out" "$sample" "${COPYKAT_CORES:-8}"
done < <("$PYTHON_BIN" - "$CKDIR" <<'PY'
import sys
import pandas as pd
ckdir = sys.argv[1]
m = pd.read_csv(f"{ckdir}/manifest.csv")
for sample in m.loc[m["status"].eq("queued"), "sample"].astype(str):
    print(sample)
PY
)
run_py s2_malignant.py --step collect
run_py s2_gse134520_breakdown.py

# ---- Stage 3: independent per-cohort clustering + MetaNeighbor raw modules ----
run_py s3_reproducibility_v2.py
run_py s3_5_characterize.py
run_py s3_6_robustness.py

# ---- Frozen signature loading + all-cell programme projection + diff_axis ----
run_py freeze_signatures.py
run_py s4_crossroute.py
run_py s4b_paired.py
if [[ "$INCLUDE_INHOUSE" == 1 ]]; then
  run_py s4c_inhouse_nodal.py
fi

# ---- Spatial + bulk portability ----
run_py s5_spatial.py
run_py s6_survival.py
run_py s6b_grade_adjust.py

# ---- Held-out validation + pooled/LOCO models + EMT diagnostic ----
run_py s12_external_validation.py
run_py s12_gse239676_validation.py
run_py s13_mixed_model.py
run_py s14_peritoneal_emt_diagnostic.py

# ---- Canonical Figure 6 source-table export ----
run_py export_fig6_jtm_source.py
