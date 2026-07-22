"""config.py — B8-GC cohort manifest + unified QC parameters + route labels.

All stages import this file. Cohort paths are resolved from environment
variables; sample-level route mappings for mixed-tissue GEO cohorts are
defined here so a clean Stage-0 rerun reproduces the locked analysis labels.

[release note] Path roots are resolved from B8GC_DATA_ROOT / B8GC_WORK_ROOT
environment variables (default: ./data, ./work relative to cwd) instead of
hardcoded absolute paths, so this config runs unmodified outside the
original deployment. Per-cohort raw-data paths below point at
{DATA_ROOT}/<cohort-specific subpath>; supply your own data under
B8GC_DATA_ROOT following the layout documented in README.md.
"""
import os
import re

HUB     = os.path.expanduser(os.environ.get("B8GC_DATA_ROOT", "./data"))
WORK    = os.path.expanduser(os.environ.get("B8GC_WORK_ROOT", "./work"))
H5AD    = f"{WORK}/h5ad"
QC      = f"{WORK}/qc"
RESULTS = f"{WORK}/results"
for d in (H5AD, QC, RESULTS):
    os.makedirs(d, exist_ok=True)

# ---- 统一 QC 阈值（全队列一致；复现性研究的前提）----
QC_PARAMS = dict(
    min_genes=500,      # 每细胞最少基因（对齐 HRA004702）
    min_counts=1000,    # 每细胞最少 UMI
    max_pct_mt=10.0,    # 线粒体比例上限
    min_cells=3,        # 每基因最少细胞
    target_sum=1e4,     # 归一化
    n_hvg=3000,         # HVG（seurat_v3, batch-aware）
    scrublet=True,      # 逐样本 doublet
)

# ---- route 分类（一等变量）----
ROUTES = ["precancer", "adjacent_normal", "primary",
          "LN_met", "ovarian_met", "peritoneal_met", "peritoneal_lavage",
          "liver_met", "ascites", "blood"]

# HRA004702 样本后缀 → route
HRA_SUFFIX2ROUTE = {
    "PT": "primary", "OM": "ovarian_met", "PM": "peritoneal_met",
    "AS": "ascites", "AD": "adjacent_normal",
}

# Verified sample-level route mappings for mixed-tissue GEO cohorts.
# GSE270680 labels "L" as lymph node without explicitly stating metastasis;
# the locked analysis uses LN_met and retains that convention here.
MAP163558 = {
    "PT": "primary", "NT": "adjacent_normal", "LN": "LN_met",
    "O": "ovarian_met", "P": "peritoneal_met", "Li": "liver_met",
}
MAP270680 = {"T": "primary", "N": "adjacent_normal", "L": "LN_met", "P": "blood"}

G183904_NORMAL = {
    "GSM5573466", "GSM5573469", "GSM5573471", "GSM5573474", "GSM5573476",
    "GSM5573486", "GSM5573488", "GSM5573490", "GSM5573496", "GSM5573500",
    "GSM5573502",  # peritoneum, histologically normal
}
G183904_PERITONEAL = {"GSM5573484", "GSM5573485", "GSM5573503"}
G183904_ALL = {f"GSM557{i}" for i in range(3466, 3506)}

G228598_LAVAGE = {
    "GSM7133742", "GSM7133743", "GSM7133744", "GSM7133747", "GSM7133749",
    "GSM7133751", "GSM7133752", "GSM7133754", "GSM7133755", "GSM7133757",
    "GSM7133758", "GSM7133759", "GSM7133760",
}
G228598_ALL = {f"GSM713{i}" for i in range(3742, 3770)}

SAMPLE_MAPPED_COHORTS = {"GSE163558", "GSE270680", "GSE183904", "GSE228598"}


def _sample_accession(sample):
    """Return the leading GSM accession from a sample/file-derived label."""
    return str(sample).split("_")[0].split(".")[0]


def route_for_sample(cohort, sample):
    """Return a verified route for a mixed-tissue cohort, or ``None`` if unknown."""
    label = str(sample)
    gsm = _sample_accession(label)
    if cohort == "GSE163558":
        match = re.match(r"([A-Za-z]+)\d*", label.split("_")[-1])
        return MAP163558.get(match.group(1) if match else None)
    if cohort == "GSE270680":
        match = re.search(r"([A-Za-z])$", label.split("_")[-1])
        return MAP270680.get(match.group(1) if match else None)
    if cohort == "GSE183904":
        if gsm in G183904_PERITONEAL:
            return "peritoneal_met"
        if gsm in G183904_NORMAL:
            return "adjacent_normal"
        if gsm in G183904_ALL:
            return "primary"
        return None
    if cohort == "GSE228598":
        if gsm in G228598_LAVAGE:
            return "peritoneal_lavage"
        if gsm in G228598_ALL:
            return "ascites"
        return None
    return None

# ---- cohort manifest ----
# fmt: mtx10x | txt | rdata | cellranger_multi
# role: discovery | validation
COHORTS = {
    "GSE163558":   dict(fmt="mtx10x", path=f"{HUB}/GC/GSE163558_raw",
                        route=None, sample_route="verified", role="discovery"),
    "GSE270680":   dict(fmt="mtx10x", path=f"{HUB}/GC/GSE270680_raw",
                        route=None, sample_route="verified", role="discovery",
                        note="EOGC; L=lymph node and is retained as LN_met to match the locked analysis"),
    "GSE228598":   dict(fmt="mtx10x", path=f"{HUB}/GC/GSE228598_raw",
                        route=None, sample_route="verified", role="discovery",
                        note="GEO sample accessions distinguish peritoneal lavage from malignant ascites"),
    "GSE134520":   dict(fmt="txt", path=f"{HUB}/GC/GSE134520_raw",
                        route="precancer", role="discovery",
                        note="NAG/CAG/IM；route 可按样本细分"),
    "GSE183904":   dict(fmt="txt", path=f"{HUB}/GC/GSE183904_raw",
                        route=None, sample_route="verified", role="discovery",
                        note="40 raw-count matrices; GEO accessions distinguish primary, adjacent-normal, "
                             "and peritoneal-tumour samples"),
    "HRA004702":   dict(fmt="rds",
                        path=f"{HUB}/HRA004702/GCscrRNAlist.rds",
                        h5ad_hint=f"{H5AD}/HRA004702.raw.h5ad",   # 转换中间文件，区别于 QDP 成品 HRA004702.h5ad
                        route=None, role="discovery",
                        sample_route="metadata",
                        route_col="orig.ident",   # 已确认：orig.ident=P01_OM/P01_PT/... ，末段后缀→route
                        note="keystone：唯一含 OM/AS。单个 Seurat(非 list)，counts 存在(~33514×281426)，"
                             "QC 状态未知 → 从 counts 重过滤；orig.ident 后缀由 _normalize_route 映射"),
    "INHOUSE_3v3": dict(fmt="rds",
                        path=f"{HUB}/GC/GC_qc.rds",
                        h5ad_hint=f"{H5AD}/INHOUSE_3v3.raw.h5ad",   # 转换中间文件，区别于 QDP 成品
                        route="primary", role="validation",
                        allow_noncounts=True,                 # 放行非整数 X（去污染分数 counts）通过 load 检查
                        preprocessed="qc_done_needs_norm",    # 已 QC(去污染/doublet)，但需按标准管线归一化
                        # 对象无 LN± 列；按 orig.ident 映射（来自你本机分析脚本的证据）
                        group_map={"P01T": "LN_pos", "P02T": "LN_pos", "P04T": "LN_pos",
                                   "P07T": "LN_neg", "P12T": "LN_neg", "P19T": "LN_neg"},
                        note="3v3 自测；单个 Seurat，6 样本。该对象从未归一化：counts 与 data 层相同，"
                             "都是去污染后的分数 counts(max≈833,非整数)，无 log 层、无原始整数 counts。"
                             "→ 导出 counts(默认)，由 s0_qdp 跳过重 QC/doublet，按标准 normalize_total+log1p 归一化"
                             "(与 discovery 一致)；route=原发，ln_status 由 group_map 注入。held-out，n=6 仅验证"),
}

DISCOVERY = [k for k, v in COHORTS.items() if v["role"] == "discovery"]
VALIDATION = [k for k, v in COHORTS.items() if v["role"] == "validation"]
