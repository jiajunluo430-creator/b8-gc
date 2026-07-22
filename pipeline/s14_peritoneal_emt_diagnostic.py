#!/usr/bin/env python
"""s14_peritoneal_emt_diagnostic.py — 诊断"HRA腹膜=非EMT vs 外部腹膜=有EMT"的矛盾。

对 HRA 内部腹膜 和 GSE308231 腹膜 做【完全相同】的 primary-vs-peritoneal 分析:
  · diff_axis Cohen's d + EMT肿瘤内在 Cohen's d (apples-to-apples,解决"是不是同一种分析")
  · 每个 EMT-TF 基因的 Cohen's d (VIM 单基因驱动? 还是 ZEB/SNAI/TWIST 真升?)
  · 每样本 EMT内在均值 (一个 outlier 样本拉的? 还是一致?)

三种可能结论:
  (a) HRA腹膜 d 也大 → Fig4c 的'flat'是欠功效/被掩盖 → 腹膜一致有EMT,Fig4c 要重看
  (b) HRA腹膜 d 小、外部大 → 真实队列异质 → 诚实写'路径依赖且队列异质'
  (c) 外部主要 VIM 升、EMT-TF 不升 → VIM 假象 → 非EMT 仍成立(VIM 广泛表达不算真EMT)
Outputs under B8GC_WORK_ROOT/qc and B8GC_WORK_ROOT/figures.
"""
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import zscore
import config as C
from fig_utils import save, W2
import s12_external_validation as s12

EMT_TF = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH2"]


def cohend(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (b.mean() - a.mean()) / sp if sp > 0 else np.nan


def prep(a, sigs):
    if "dax" not in a.obs:
        for k in ["M1", "M2"]:
            g = [x for x in sigs[k] if x in a.var_names]
            if g:
                sc.tl.score_genes(a, g, score_name=f"s_{k}")
        a.obs["dax"] = zscore(a.obs["s_M1"].values) - zscore(a.obs["s_M2"].values)
    egi = [g for g in EMT_TF if g in a.var_names]
    sc.tl.score_genes(a, egi, score_name="s_EMTi")
    return egi


def analyze(a, label, sigs):
    egi = prep(a, sigs)
    rt = a.obs["route"].astype(str)
    pmask = rt.str.contains("primary|prim|tumor|in_situ|insitu", case=False, na=False).values
    qmask = rt.str.contains("periton", case=False, na=False).values
    npri, nper = int(pmask.sum()), int(qmask.sum())
    if npri < 20 or nper < 20:
        return dict(label=label, ok=False, npri=npri, nper=nper)
    dax = a.obs["dax"].values; emti = a.obs["s_EMTi"].values
    d_diff = cohend(dax[pmask], dax[qmask]); d_emti = cohend(emti[pmask], emti[qmask])
    sub = a[:, egi]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    gd = {g: cohend(X[pmask, i], X[qmask, i]) for i, g in enumerate(egi)}
    scol = next((c for c in a.obs.columns if c.lower() in
                 ("sample", "patient", "orig.ident", "sample_id", "patient_id", "samplename")), None)
    sp = sm = None
    if scol:
        s = pd.DataFrame({"s": a.obs[scol].astype(str).values, "e": emti, "p": pmask, "q": qmask})
        sp = s[s.p].groupby("s")["e"].mean(); sm = s[s.q].groupby("s")["e"].mean()
    return dict(label=label, ok=True, npri=npri, nper=nper, d_diff=d_diff, d_emti=d_emti, gd=gd, sp=sp, sm=sm)


def main():
    sigs = json.load(open(f"{C.RESULTS}/signatures.json"))
    L = ["# S14 腹膜 EMT 矛盾诊断(HRA内部腹膜 vs GSE308231腹膜,同口径)\n"]
    res = []

    hpath = f"{C.RESULTS}/S4_scored.h5ad"
    if os.path.exists(hpath):
        h = sc.read_h5ad(hpath)
        L.append(f"- HRA: {hpath} 载入 {h.n_obs} 细胞; route 值={sorted(set(h.obs['route'].astype(str)))[:8]}")
        res.append(analyze(h, "HRA peritoneal", sigs))
    else:
        L.append(f"⚠ 未找到 HRA {hpath} — 贴我 S4 输出的 h5ad 路径")

    cfg = next(c for c in s12.COHORTS if c["name"] == "GSE308231")
    m, diag, files = s12.load_cohort(cfg)
    if m is not None:
        e, cl, keep = s12.epithelial(m)
        res.append(analyze(e, "GSE308231 peritoneal", sigs))
    else:
        L.append("⚠ GSE308231 载入失败: " + "; ".join(diag))

    for r in res:
        if not r.get("ok"):
            L.append(f"## {r['label']}: 样本不足(primary {r.get('npri')} / peri {r.get('nper')} 细胞)\n"); continue
        L.append(f"## {r['label']}  (primary {r['npri']} / peritoneal {r['nper']} 细胞)")
        L.append(f"- **diff_axis d={r['d_diff']:+.2f} | EMT肿瘤内在 d={r['d_emti']:+.2f}**")
        L.append("- 每基因 d(peri vs primary): " + ", ".join(f"{g}={d:+.2f}" for g, d in r["gd"].items()))
        if r["sp"] is not None:
            L.append(f"- 每样本 EMT内在: primary={list(np.round(r['sp'].values,2))} vs peri={list(np.round(r['sm'].values,2))}")
        L.append("")

    if len(res) == 2 and all(r.get("ok") for r in res):
        h, g = res[0], res[1]
        L.append("## 对比结论")
        L.append(f"- EMT肿瘤内在 d: HRA腹膜={h['d_emti']:+.2f} vs GSE308231腹膜={g['d_emti']:+.2f}")
        if abs(h["d_emti"]) < 0.4 and abs(g["d_emti"]) >= 0.5:
            L.append("→ (b) **真实队列异质**: HRA腹膜非EMT、外部腹膜有EMT,不是 bug;诚实写'EMT 参与路径依赖且队列异质'。")
        elif abs(h["d_emti"]) >= 0.5:
            L.append("→ (a) **HRA腹膜其实也EMT**: Fig4c 的'flat'是欠功效/被均值掩盖 → 腹膜一致有EMT,Fig4c 那条线要重看。")
        else:
            L.append("→ 两边都不大 → 腹膜 EMT 信号弱,非EMT 大体可守。")
        tfs = ["ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1"]
        gtf = np.nanmean([abs(g["gd"].get(t, np.nan)) for t in tfs])
        gvim = abs(g["gd"].get("VIM", np.nan))
        L.append(f"- GSE308231 驱动基因: VIM d={g['gd'].get('VIM', float('nan')):+.2f} vs EMT-TF 均|d|={gtf:.2f} → "
                 + ("(c 排除) EMT-TF 真升 = genuine EMT" if gtf > 0.3 else "(c) 主要 VIM 驱动 → 广泛表达,真 EMT 证据弱,非EMT 仍可守"))

    oks = [r for r in res if r.get("ok")]
    if oks:
        genes = [g for g in EMT_TF if all(g in r["gd"] for r in oks)]
        fig, ax = plt.subplots(figsize=(W2, 0.4 * W2))
        x = np.arange(len(genes)); wd = 0.8 / len(oks)
        for j, r in enumerate(oks):
            ax.bar(x + j * wd, [r["gd"][g] for g in genes], wd, label=r["label"], color=["#457b9d", "#e76f51"][j % 2])
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_xticks(x + wd * (len(oks) - 1) / 2); ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=6)
        ax.set_ylabel("Cohen's d (peritoneal vs primary)"); ax.legend(fontsize=6)
        ax.set_title("Per-gene EMT change: HRA vs external peritoneal")
        plt.tight_layout(); save(fig, "FigS8_peri_emt")

    open(f"{C.QC}/S14_peri_emt.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
