#!/usr/bin/env python
"""s2_gse134520_breakdown.py — sample/tissue-level malignant breakdown for GSE134520."""
import re
import pandas as pd
import scanpy as sc
import config as C


def tissue_class(sample):
    s = str(sample).upper()
    for k in ["NAG", "CAG", "IMW", "IMS", "EGC"]:
        if k in s:
            return k
    return "OTHER"


def find_ctcol(obs):
    for c in obs.columns:
        vals = " ".join(str(x).lower() for x in pd.unique(obs[c])[:40])
        if "epitheli" in vals and ("myeloid" in vals or "fibro" in vals or "t_nk" in vals or "plasma" in vals):
            return c
    return "cell_type" if "cell_type" in obs.columns else obs.columns[0]


def main():
    s1 = sc.read_h5ad(f"{C.RESULTS}/S1_integrated.h5ad", backed="r").obs.copy()
    s2 = sc.read_h5ad(f"{C.RESULTS}/S2_malignant.h5ad", backed="r").obs.copy()
    s1 = s1[s1["cohort"].astype(str) == "GSE134520"].copy()
    s2 = s2[s2["cohort"].astype(str) == "GSE134520"].copy()
    ctcol = find_ctcol(s1)
    epi = s1[s1[ctcol].astype(str).eq("Epithelial")].copy()
    epi["tissue_class"] = epi["sample"].map(tissue_class)
    s2["tissue_class"] = s2["sample"].map(tissue_class)

    epi_n = epi.groupby("sample", observed=True).size().rename("n_epithelial")
    mal_n = s2.groupby("sample", observed=True).size().rename("n_malignant")
    tab = pd.concat([epi_n, mal_n], axis=1).fillna(0).astype(int).reset_index()
    tab["tissue_class"] = tab["sample"].map(tissue_class)
    tab["malignant_fraction"] = tab["n_malignant"] / tab["n_epithelial"].replace(0, pd.NA)
    tab = tab.sort_values(["tissue_class", "sample"])

    agg = (tab.groupby("tissue_class", observed=True)[["n_epithelial", "n_malignant"]]
             .sum()
             .reset_index())
    agg["malignant_fraction"] = agg["n_malignant"] / agg["n_epithelial"].replace(0, pd.NA)

    tab.to_csv(f"{C.RESULTS}/S2_GSE134520_breakdown_by_sample.csv", index=False)
    agg.to_csv(f"{C.RESULTS}/S2_GSE134520_breakdown_by_tissue.csv", index=False)

    L = ["# S2 GSE134520 malignant breakdown",
         f"- epithelial cells: {epi.shape[0]}",
         f"- malignant calls: {s2.shape[0]}",
         f"- cell_type column used: {ctcol}",
         "",
         "## by tissue class",
         agg.to_markdown(index=False),
         "",
         "## by sample",
         tab.to_markdown(index=False)]
    with open(f"{C.QC}/S2_GSE134520_breakdown.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
