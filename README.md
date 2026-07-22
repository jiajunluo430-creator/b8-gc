# B8-GC: cross-cohort reproducible malignant cell-state atlas of gastric cancer metastatic routes

This repository is a public, sanitised code release accompanying the B8-GC
manuscript. It provides the code needed to reproduce the manuscript's core methodological stages —
QC/atlas construction, CopyKAT-based malignant identification, independent
per-cohort clustering with MetaNeighbor cross-cohort validation, frozen
signature projection, differentiation-axis scoring, route/patient-paired
summaries, EMT scoring, spatial analysis, held-out validation, and bulk
survival portability — from permitted, non-identifiable inputs.

**Locked submission commit** (of the internal analysis repository this
release was exported from): `fcbf48d94a8dfd806d6d4a085f60289c0dcb012f`. See
`docs/PROVENANCE.md` for full provenance, frozen values, and prohibited
claims.

## Study overview

Malignant epithelial cell states across independent gastric cancer scRNA-seq
cohorts are identified by copy-number aneuploidy (CopyKAT), clustered
independently per cohort, and tested for cross-cohort reproducibility with
MetaNeighbor. Three frozen 30-gene programme signatures (M0/M1/M2) — derived
once from this raw-module structure and never re-fit — are then projected
onto *every* carcinoma-derived malignant cell, regardless of which raw module
(or none) that cell's cluster belonged to. The resulting M1–M2 differentiation
axis is characterised across metastatic routes, patient-paired primary/met
samples, spatial transcriptomics, held-out external cohorts, and bulk
survival data.

### The two-stage raw-module -> frozen-programme framework

This is the single most important structural idea in the codebase, and the
one most likely to be misread if skipped:

1. **Stage 1 — raw module discovery** (`pipeline/s3_reproducibility_v2.py`,
   `pipeline/run_metaneighbor.R`, `pipeline/s3_5_characterize.py`,
   `pipeline/s3_6_robustness.py`). Clusters are discovered **independently
   within each cohort** (no cross-cohort integration at this step), then
   tested for cross-cohort replication with MetaNeighbor at a locked graph
   threshold (edge AUROC = 0.86, requiring recovery in >=3 cohorts). This
   currently yields **five** raw reproducible cross-cohort modules,
   **RMC0–RMC4**, plus a residual `unrepro` set of clusters that do not reach
   the 3-cohort criterion.

2. **Stage 2 — frozen programme projection** (`pipeline/sig_utils.py`,
   `pipeline/freeze_signatures.py`, `pipeline/s4_crossroute.py`). Frozen
   30-gene one-versus-rest differential-expression signatures for **three**
   programmes — **M0, M1, M2** — were derived once from three of the five raw
   modules and are applied unchanged (argmax of per-cell signature scores) to
   **every** carcinoma-derived malignant cell, including cells from RMC0,
   RMC2, or `unrepro`. No cell is excluded from Stage-2 scoring on the basis
   of its Stage-1 module membership.

**Current five RMCs vs. final three programmes — do not conflate these two
token spaces:**

| Raw module (Stage 1) | Frozen programme (Stage 2) |
|---|---|
| RMC1 | -> M0 |
| RMC3 | -> M1 |
| RMC4 | -> M2 |
| RMC0, RMC2, `unrepro` | *(no dedicated programme — still scored/argmax'd into M0/M1/M2)* |

The differentiation axis reported throughout the manuscript is
`diff_axis = s_M1 - s_M2` (`pipeline/s4_crossroute.py`); M0 is an auxiliary
pit/metabolic facet, not a pole of this axis.

See `docs/PROVENANCE.md` for the full list of permitted/prohibited claims
around this framework (in particular: **do not** report the Stage-1 AUROC
threshold (0.86) as a Stage-2/final-programme AUROC).

## Installation

```bash
# Exact Linux x86-64 analysis environment:
conda env create -f environment.yml
conda activate b8gc

# Or the complete pip environment in an existing Python 3.11 installation:
pip install -r requirements.txt

# Lightweight smoke-test environment:
pip install -r requirements-test.txt
python -m pytest -q
```

`environment.yml` is the exact **linux-64** lock captured for the submitted
analysis; its build strings are not intended as a cross-platform lock.
`requirements-test.txt` contains only the Python packages needed for the
included tests and synthetic demonstration. GitHub Actions repeats the Python,
Bash, and R syntax/smoke checks on every push and pull request.

R stages (`export_rds.R`, `run_copykat.R`, `run_metaneighbor.R`) require R
with `Seurat`, `Matrix`, `copykat`, and `MetaNeighbor` (Bioconductor)
installed. See `docs/environment/` for the exact package versions used to
generate the manuscript's results.

## Expected input formats

Raw cohort data are **not included** in this repository (see "Data access"
below). Point `B8GC_DATA_ROOT` (default `./data`) at a directory laid out as:

```
$B8GC_DATA_ROOT/
├── GC/GSE163558_raw/            # 10x mtx triplets, one dir per sample
├── GC/GSE270680_raw/            # 10x mtx triplets
├── GC/GSE228598_raw/            # 10x mtx triplets
├── GC/GSE134520_raw/            # tab-delimited count matrices
├── GC/GSE183904_raw/            # per-sample *.csv.gz (genes x cells)
├── HRA004702/GCscrRNAlist.rds   # Seurat object, GSA-accessioned
├── GC/GC_qc.rds                 # optional in-house 3-vs-3 cohort, not bundled
├── GC/TCGA_STAD_Exp.txt, TCGA_STAD_time.txt, TCGA_STAD_clinical.txt
├── GC/GSE84437_raw/GSE84437_estimate_input.txt, GSE84437_survival_raw.csv
├── GC/GSE251950_raw/             # Visium: *_filtered_feature_bc_matrix.h5 + tissue_positions
├── GSE239676_raw/                # 10x mtx.gz + barcodes/features/meta.tsv.gz
└── external/
    ├── GSE308231/GSE308231_RAW.tar
    └── GSE246662/GSE246662_RAW.tar
```

Exact per-cohort paths and formats are declared in `pipeline/config.py`
(`COHORTS` dict) and the per-script constants in `pipeline/s6_survival.py`,
`pipeline/s6b_grade_adjust.py`, `pipeline/s5_spatial.py`,
`pipeline/s12_external_validation.py`, and `pipeline/s12_gse239676_validation.py`.
All of these resolve through `B8GC_DATA_ROOT` (or dedicated
`B8GC_<COHORT>_ROOT` overrides) rather than hardcoded absolute paths.
The four mixed-tissue discovery cohorts are assigned at sample level by the
verified mappings in `pipeline/config.py`; unexpected sample labels are marked
`unknown` instead of silently receiving a cohort-wide route.

## Stage-by-stage workflow

| Stage | Scripts | Command |
|---|---|---|
| Unified QC + atlas construction | `config.py`, `s0_qdp.py`, `s1_integrate.py`, `build_h5ad.py`, `export_rds.R` | `python pipeline/s0_qdp.py --cohort <name>`, then `python pipeline/s1_integrate.py` |
| Per-sample CopyKAT malignant ID | `s2_malignant.py`, `s2_gse134520_breakdown.py`, `run_copykat.R` | `python pipeline/s2_malignant.py --step prep`, run CopyKAT per sample, `--step collect` |
| Independent clustering + MetaNeighbor raw modules | `s3_reproducibility_v2.py`, `run_metaneighbor.R`, `s3_5_characterize.py`, `s3_6_robustness.py` | `python pipeline/s3_reproducibility_v2.py` |
| Frozen signature loading + all-cell programme projection + differentiation-axis scoring | `sig_utils.py`, `freeze_signatures.py`, `s4_crossroute.py` | `python pipeline/freeze_signatures.py && python pipeline/s4_crossroute.py` |
| Route/sample/patient-paired summaries | `s4b_paired.py`, `s4c_inhouse_nodal.py` | `python pipeline/s4b_paired.py` |
| EMT scoring | `s14_peritoneal_emt_diagnostic.py` | `python pipeline/s14_peritoneal_emt_diagnostic.py` |
| Spatial Moran analysis | `s5_spatial.py` | `python pipeline/s5_spatial.py` |
| Held-out validation + pooled/LOCO models | `s12_external_validation.py`, `s12_gse239676_validation.py`, `s13_mixed_model.py` | `python pipeline/s12_external_validation.py` |
| Bulk negative-portability survival | `s6_survival.py`, `s6b_grade_adjust.py` | `python pipeline/s6_survival.py` |
| Figure/source-table export | `fig_utils.py`, `style.py`, `export_fig6_jtm_source.py` | `python pipeline/export_fig6_jtm_source.py` |

Or run the whole ordered sequence with `./run_all.sh` (set `B8GC_DATA_ROOT` /
`B8GC_WORK_ROOT` first, or rely on the `./data` / `./work` defaults). The
in-house stage is included automatically only when `GC/GC_qc.rds` is present;
the public-cohort workflow otherwise proceeds without it.

For a self-contained demonstration of the Stage-2 frozen-signature projection
and differentiation-axis logic on fully synthetic data (no cohort data
required), see `examples/synthetic_demo/`.

## Scope of this release

The 11 categories above are a **methodologically-focused subset** of the full
internal analysis pipeline. The full pipeline additionally includes SCENIC/TF
activity, gene-knockout modeling, trajectory analysis, pathway enrichment,
Lauren-classification panels, supplementary CNV analyses, and the full
multi-panel figure-rendering suite — these are exploratory/stylistic and are
out of scope for this release.

Three scripts were deliberately **excluded** as retired/non-canonical, not
silently dropped:

- **`s3_reproducibility.py` (v1)** — an earlier, methodologically retired
  approach that defined malignant states *inside* an already
  Harmony-integrated embedding and then tested cross-cohort AUROC on that
  same integrated space — a circular test, since integration already forces
  cross-cohort alignment (cross-cohort AUROC trivially approaches 1). It also
  mixes in technical artifacts (e.g. stress/OxPhos response) as if they were
  biological states. Its successor, `s3_reproducibility_v2.py` (independent
  per-cohort clustering, *then* cross-cohort MetaNeighbor testing), is the
  only version any downstream script reads from.
- **`remap_routes.py`** — an uncalled, standalone patch script that is not
  invoked by `run_all.sh` or any other script in the pipeline; its
  verified route-mapping logic is encoded directly in `pipeline/config.py` /
  `pipeline/s0_qdp.py` and covered by `tests/test_route_mapping.py`.
- **`emit_manuscript_values.py`** — an internal manuscript-assembly helper
  that mixed live results with hardcoded headline counts and out-of-scope
  rendering/trajectory artifacts. It cannot be regenerated solely from this
  release and is therefore not represented as part of the public pipeline.
  Canonical Figure 6 source tables are generated directly from the Stage 12/13
  outputs by `export_fig6_jtm_source.py`.

Also excluded: the retired Figure 2 "three-state plateau"
(`Fig2_C_threshold_sweep.csv`, a hardcoded literal array disconnected from
any real computation, and the two render scripts that consumed it) — see
`docs/PROVENANCE.md` for the full account of why this panel was retired.

## Data access

- **Public records/accessions**: GSE163558, GSE270680, GSE228598, GSE134520,
  GSE183904, GSE251950, GSE239676, GSE308231, GSE246662, GSE84437, and
  TCGA-STAD are deposited in GEO / GDC. Raw-versus-processed availability
  differs by accession; follow each repository record and arrange the files in
  the layout above. In particular, some human cohorts provide processed
  matrices while restricting raw sequence data.
- **HRA004702** is an open-access GSA-Human accession.
- **`INHOUSE_3v3`** (the in-house 3-vs-3 lymph-node-status validation cohort)
  is not included in this repository. Its availability is governed by the
  final Data Availability statement of the associated manuscript.

## Citation

If you use this code, cite the associated manuscript:

> Luo J, Wang Z, Liang X. *A reproducible gastric differentiation axis marks
> metastatic dedifferentiation in gastric cancer.* 2026.

Repository: https://github.com/jiajunluo430-creator/b8-gc

Archived release: https://doi.org/10.5281/zenodo.21496971

The root `CITATION.cff` supplies machine-readable citation metadata to GitHub.

## Repository layout

```
b8-gc/
├── README.md                 # this file
├── LICENSE                   # MIT
├── CITATION.cff              # repository/manuscript citation metadata
├── environment.yml           # conda env (Python 3.11.15 + pinned pip deps)
├── requirements.txt          # complete locked pip environment
├── requirements-test.txt     # lightweight smoke-test environment
├── run_all.sh                # canonical ordered re-run entrypoint
├── docs/
│   ├── PROVENANCE.md         # lock commit, RMC/programme crosswalk, frozen values, prohibited claims
│   └── environment/          # R session info, curated R package list, pip freeze
├── signatures/                # frozen M0/M1/M2 gene lists + checksums (see signatures/SHA256SUMS.txt)
├── pipeline/                  # the sanitised analysis scripts (see table above)
├── examples/synthetic_demo/   # synthetic end-to-end demo, no real data required
├── tests/                     # pytest smoke tests
├── RELEASE_MANIFEST.tsv
├── QC_REPORT.md
└── SHA256SUMS.txt
```
