# Provenance

## Submission lock

This release was exported from the internal B8-GC analysis repository at
commit `fcbf48d94a8dfd806d6d4a085f60289c0dcb012f` (the analysis lock used for
the pre-submission manuscript package). No tracked file in the internal repository was
modified to produce this release; the release is an export/sanitisation of
that locked state into a new directory.

## The raw-module / frozen-programme distinction

This is the terminology contract for the whole codebase and manuscript.
Getting it wrong inverts what the reproducibility analysis actually claims.

### Stage 1 — raw reproducible cross-cohort modules (RMCs)

Computed by `pipeline/s3_reproducibility_v2.py`:

1. Malignant (CopyKAT-called aneuploid) cells are clustered **independently
   within each cohort** — no Harmony/cross-cohort integration is applied
   before this clustering step.
2. Per-cohort clusters are cross-checked pairwise with **MetaNeighbor**
   (`pipeline/run_metaneighbor.R`), producing a cluster x cluster AUROC
   matrix.
3. An edge is drawn between two clusters (from different cohorts) when their
   MetaNeighbor AUROC >= **0.86**.
4. A connected component of mutually-replicating clusters is called a
   **reproducible module** if it is recovered in **>= 3 cohorts**.

At the locked commit, this procedure currently yields **five** such modules:
**RMC0, RMC1, RMC2, RMC3, RMC4**, plus a residual `unrepro` bucket for
clusters that fail the >=3-cohort criterion. RMC indices are stable
identifiers used throughout `results/` and QC summaries — they are **not**
biological programme names.

**Prohibited claim**: do not describe the 0.86 MetaNeighbor edge threshold,
or any RMC-level AUROC, as the reproducibility statistic of the *final*
M0/M1/M2 programmes. The edge threshold is a Stage-1 graph-construction
parameter over raw per-cohort clusters, not a property of the frozen
signatures.

### Stage 2 — frozen programme signatures (M0/M1/M2)

Computed by `pipeline/freeze_signatures.py` (signature derivation, one-time)
and applied by `pipeline/sig_utils.py` / `pipeline/s4_crossroute.py`
(signature scoring, applied unchanged thereafter):

- Three 30-gene one-versus-rest differential-expression signatures — **M0**,
  **M1**, **M2** — were derived once, from three of the five Stage-1 raw
  modules, and then frozen: `signatures/signatures.json` /
  `signatures/SuppTable_S2_signatures.csv` are the exact, final gene lists,
  locked by the checksum in `signatures/SHA256SUMS.txt` and enforced at
  runtime via `assert_signature_hash()` / `EXPECTED_SIGNATURE_SHA256` in
  `pipeline/sig_utils.py`.
- The frozen signatures are then projected, by `argmax` of per-cell
  signature score, onto **every** carcinoma-derived malignant cell in every
  cohort — including cells whose Stage-1 cluster fell in RMC0, RMC2, or
  `unrepro`. Stage-1 module membership is not a gate on Stage-2 scoring.

### Crosswalk (Stage 1 raw module -> Stage 2 frozen programme)

| Raw module (RMC) | Frozen programme |
|---|---|
| RMC1 | M0 |
| RMC3 | M1 |
| RMC4 | M2 |
| RMC0, RMC2, `unrepro` | *(not a signature source; still scored into M0/M1/M2)* |

`diff_axis = s_M1 - s_M2` (`pipeline/s4_crossroute.py`) is the
differentiation axis reported throughout the manuscript's route/paired/
spatial/validation/survival analyses. M0 is a pit/metabolic auxiliary facet
and is not one of the two poles of `diff_axis`.

## Route-label source records

Sample-level route assignments for the four mixed-tissue discovery cohorts
are locked in `pipeline/config.py` and trace to their official GEO records:

- [GSE163558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163558):
  primary, adjacent normal, lymph node, ovary, peritoneum, and liver samples.
- [GSE270680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270680):
  tumour, normal, lymph-node, and peripheral-blood samples. GEO labels the
  `L` samples as lymph node without explicitly stating metastatic status; the
  locked analysis retains `LN_met` and this caveat must accompany that label.
- [GSE183904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183904):
  primary gastric tumour/normal and peritoneal tumour/normal samples, assigned
  by GSM accession.
- [GSE228598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228598):
  malignant ascites and peritoneal-washing samples, assigned by GSM accession.

Unexpected sample labels return `unknown`; they are never silently assigned a
cohort-wide route. `tests/test_route_mapping.py` locks representative mappings
and the unknown-label guard.

## Retired / excluded analyses

### Retired: `s3_reproducibility.py` (v1)

An earlier approach that defined malignant cell states **inside an
already-Harmony-integrated embedding**, then computed cross-cohort
MetaNeighbor AUROC on that same integrated space. This is circular: Harmony
integration is explicitly designed to force cross-cohort alignment, so
testing "cross-cohort reproducibility" on states defined post-integration
trivially drives AUROC toward 1 regardless of whether the underlying
biology replicates. Several of its resulting "states" are also technical
(stress-response / OxPhos) rather than biological. This script, and its
output `S3_malignant_states.h5ad`, are read by nothing downstream and are
**excluded** from this release. `s3_reproducibility_v2.py` — independent
per-cohort clustering *before* any cross-cohort comparison — is the only
version any other script in this release reads from, and is the canonical
Stage-1 analysis.

### Retired: Figure 2 "three-state plateau"

The original three-malignant-state plateau figure was driven by
`figure_source_data/Fig2_C_threshold_sweep.csv`, a hardcoded literal array
disconnected from any live computation in the pipeline, rendered by
`Fig2_render.R` / `b8_pipeline/render_batch1_fig1_2_6.{py,R}`. None of these
three files are included in this release. The current Figure 2 (see the
canonical export process referenced in the manuscript) supersedes this
retired panel.

### Excluded as evidence: Phase 1B

Any "Phase 1B" intermediate audit materials referenced in internal working
notes are preliminary/exploratory checkpoints, not submission evidence, and
are not part of this release.

### Excluded after integration: `remap_routes.py`

A standalone patch script whose docstring references a private
manual-verification document. The sample-level mappings it established for
GSE163558, GSE270680, GSE183904, and GSE228598 are integrated into
`pipeline/config.py` / `pipeline/s0_qdp.py` and covered by
`tests/test_route_mapping.py`; the one-off mutating patch itself is excluded.

### Excluded as non-self-contained: `emit_manuscript_values.py`

The internal manuscript-assembly helper mixed live pipeline outputs with
hardcoded cohort counts and figure/trajectory/quarantine artifacts outside the
scope of this release. Including it would make `run_all.sh` fail even after all
documented inputs were supplied. The helper is therefore excluded rather than
presented as reproducible public code. Figure 6 source tables are instead
derived directly from Stage 12 cell tables and the Stage 13 pooled-model JSON
by `pipeline/export_fig6_jtm_source.py`.

## Primary frozen values

The headline pooled differentiation-axis effect reported in Figure 6 is
computed by `pipeline/s13_mixed_model.py` and written to
`results/S13_pooled.json` -> `validation_only`
(the three-cohort held-out validation-only pooled model: GSE308231 (peritoneal),
GSE246662 (liver), GSE239676 (liver/ascites)):

```
beta = -1.0097119609277696
CI   = [-1.498143203859601, -0.5212807179959382]
P    = 0.0001396455390576807
```

This is a *validation-only* pooled model (independent held-out cohorts),
distinct from the `all_cohorts`, `validation_fixed_meta`, and `fixed_meta`
values also present in `S13_pooled.json`, some of which pool across cohorts used
in earlier stages of the analysis. When citing "the" Figure 6 pooled effect,
use the value above unless the manuscript text specifies otherwise.

## Bulk negative-portability result

`pipeline/s6_survival.py` / `pipeline/s6b_grade_adjust.py` test whether
`diff_axis`, scored on bulk RNA-seq deconvolution proxies, portrays any
association with overall survival in two independent bulk cohorts:

- **TCGA-STAD** (n=352): Cox HR=0.967, p=0.7753; KM log-rank p=0.4093.
- **GSE84437** (n=431): Cox HR=1.194, p=0.2808; KM log-rank p=0.0192.

Both are reported as **no significant association** — this is a deliberate,
honestly-reported negative result establishing that `diff_axis` does not
naively transfer to bulk deconvolution proxies, not a result to be
re-interpreted as a positive finding by post-hoc subgroup search.
