# QC report — B8-GC public code release

Audit date: 2026-07-22
Internal analysis lock reported by the source package:
`fcbf48d94a8dfd806d6d4a085f60289c0dcb012f`

## 1. Audit scope and baseline integrity

The supplied archive was reviewed in an isolated extraction before any edits.
Its SHA256 was:

```
036ab861d8c95f81a96a0b7dc92ac4ea44494d3df0984e5d538e9c071f1803ab
```

Baseline archive checks:

- 53 files under one top-level release directory; no absolute, parent-traversal,
  duplicate, or symbolic-link entries.
- All 52 entries in the supplied root `SHA256SUMS.txt` matched; the only file
  not listed was the checksum file itself.
- Both frozen-signature file hashes matched `signatures/SHA256SUMS.txt`.
- The canonical semantic signature hash matched the runtime lock:
  `b522d543b269a86476499ee2955c35a90cdcd3bce68c9d3f5f9c3634f7377324`.

The supplied archive was not modified. All corrections were made in a separate
GitHub-ready working copy. The internal server repository state was not
independently re-queried during this Windows-side audit.

## 2. Material corrections made during independent review

- Integrated verified sample-level route assignments for GSE163558,
  GSE270680, GSE183904, and GSE228598 into the canonical Stage-0 path. The
  previous cohort-wide defaults would have mislabeled mixed-tissue samples on
  a clean rerun. Unknown sample labels now fail safe to `unknown`.
- Implemented actual leave-one-validation-cohort-out fixed-effect analyses in
  `s13_mixed_model.py`, while keeping the headline validation-only
  cohort-adjusted sample-level OLS distinct from fixed-effect meta-analysis.
- Rebuilt `export_fig6_jtm_source.py` to read only outputs created by the
  released Stage 12/13 scripts. Six unpublished legacy Figure 6 intermediates
  are no longer required.
- Removed `emit_manuscript_values.py`: it mixed live results with hardcoded
  counts and out-of-scope internal artifacts, so it could not be regenerated
  from the public release.
- Made the in-house 3-vs-3 stage optional in `run_all.sh`; a public-cohort run
  no longer fails solely because the non-bundled in-house RDS is absent.
- Corrected the signature checksum verification command and the malformed
  `packaging==26.2` entry in the environment snapshot.
- Added a minimal test environment, GitHub Actions checks, route/LOCO/exporter
  regression tests, and schema-valid citation metadata.

## 3. Privacy and repository-hygiene sweep

A full-tree sweep after the corrections found:

- no author-specific usernames or user-home paths;
- no unresolved bracketed author/repository/data-access placeholders;
- no credential-shaped assignments, private-key blocks, or obvious secrets;
- no raw or processed patient/cohort matrices, sequencing files, or generated
  analysis objects;
- no cache directories, bytecode, test outputs, or synthetic-demo outputs.

The final tree contains 58 small code/documentation files (approximately
251 KB before the checksum file is regenerated); the largest file is a Python
script, not a data artifact. All cohort and scratch paths use `B8GC_*`
environment variables with relative defaults.

## 4. Independent verification results

Static and metadata checks:

- 24 `pipeline/*.py` scripts and all tests/examples compile cleanly.
- All 3 `pipeline/*.R` scripts parse under R 4.5.2.
- `run_all.sh` passes `bash -n`.
- `CITATION.cff` validates against Citation File Format schema 1.2.0.
- Frozen JSON and supplementary-table signatures still round-trip to the same
  three 30-gene programmes and canonical semantic hash.

Python smoke/regression suite, run in a disposable environment created from
`requirements-test.txt`:

```
21 passed, 115 warnings in 29.01s
```

The warnings were dependency deprecation and AnnData index-conversion warnings;
there were no failed assertions. Added tests cover:

- all-cell Stage-2 denominator retention;
- signature counts, uniqueness, tamper detection, and canonical hash;
- supplementary/source-table schemas;
- sample-level GEO route mappings and unknown-label guards;
- locked fixed-effect and LOCO estimates;
- Figure 6 source export from live Stage 12/13 schemas and absence of legacy
  Figure 6 input dependencies.

Locked numerical regression results:

| Analysis | Estimate | 95% CI |
|---|---:|---:|
| Validation-only sample-level OLS (reference) | -1.0097119609 | -1.4981432039 to -0.5212807180 |
| Validation fixed-effect meta-analysis | -1.2316571205 | -1.6558120797 to -0.8075021613 |
| LOCO: omit GSE308231 | -1.2870233105 | -1.7338580885 to -0.8401885325 |
| LOCO: omit GSE246662 | -0.8820827423 | -1.4041564473 to -0.3600090373 |
| LOCO: omit GSE239676 | -1.6437407081 | -2.2839982178 to -1.0034831985 |

The fully synthetic 2,880-cell by 208-gene demonstration also ran end to end:

- primary mean `diff_axis` 0.254;
- ovarian metastasis -0.039;
- peritoneal metastasis -0.121;
- ascites -0.239;
- simulated-programme argmax agreement 90.8% (intentionally imperfect).

## 5. Verification boundary

No real cohort data are bundled. Consequently, Stage 0 cohort QC, CopyKAT,
MetaNeighbor, spatial, survival, and external-validation stages were not rerun
end to end in this audit. Syntax checks and unit/synthetic tests do not prove
that every data-dependent branch executes on every source dataset. The report
therefore makes no claim that all 11 categories were freshly reproduced.

The four corrected route mappings were checked against official GEO sample
metadata and are documented in `docs/PROVENANCE.md`. GSE270680's `L` samples
are described by GEO as lymph node without explicit metastatic wording; the
locked `LN_met` convention is retained with that caveat.

## 6. Public-release metadata

- `CITATION.cff` uses the manuscript's ordered authors and affiliations:
  Jiajun Luo, Ziwei Wang, and Xiaolong Liang.
- The MIT notice names the same three authors as copyright holders.
- The repository URL is `https://github.com/jiajunluo430-creator/b8-gc`.
- The versioned Zenodo archive DOI is
  `https://doi.org/10.5281/zenodo.21496971`.
- The manuscript's final Data Availability statement remains the authority for
  the in-house cohort and must remain consistent with the neutral README wording.
