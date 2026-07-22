# Synthetic demo: Stage-2 frozen-signature projection

Demonstrates the Stage-2 logic from `pipeline/s4_crossroute.py` — frozen
M0/M1/M2 signature scoring, argmax programme assignment, and
`diff_axis = s_M1 - s_M2` — on a small, fully synthetic dataset. No real or
patient data is used anywhere in this demo.

`generate_synthetic_data.py` builds an AnnData over the 90 frozen signature
genes plus 120 filler genes, with 3 fake cohorts x 4 fake routes x 4 fake
samples x 60 cells, where each cell's latent programme is drawn from a
route-dependent distribution (metastatic routes skewed toward M2) purely so
the demo has a directional pattern worth recovering.

## Run

```bash
cd examples/synthetic_demo
python generate_synthetic_data.py
python run_demo.py
```

Expected output: a per-route argmax programme composition table, a
per-route mean `diff_axis` table (primary should show the highest mean
`diff_axis`, `ascites` the lowest, reflecting the simulated skew), an
argmax-vs-ground-truth agreement rate (not expected to be 100% — signature
scoring is inherently noisy), and `synthetic_demo_summary.csv`.
