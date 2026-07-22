#!/usr/bin/env Rscript
# 用法: Rscript export_rds.R <input.rds> <outdir> [layer]
#   layer 可选，默认 "counts"；若原始 counts 不可恢复(被去污染/SCT覆盖)，可传 "data"(log归一化)。
# 仅依赖 Seurat + Matrix。导出所选层(genes×cells) + meta → build_h5ad.py 组装 h5ad(放入 X)。
suppressMessages({library(Seurat); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 2)
obj <- readRDS(a[1]); outdir <- a[2]
layer <- if (length(a) >= 3) a[3] else "counts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (is.list(obj) && !inherits(obj, "Seurat")) obj <- merge(obj[[1]], y = obj[-1])
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)   # Seurat v5 多层→合并；v4 无操作

# --- 诊断：对象里到底有什么 ---
cat("assays:", paste(tryCatch(SeuratObject::Assays(obj), error=function(e) names(obj@assays)),
                     collapse=","), "\n")
cat("RNA layers:", paste(tryCatch(SeuratObject::Layers(obj[["RNA"]]),
                                  error=function(e) "(v4: counts/data/scale.data)"), collapse=","), "\n")

mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = layer),
                error = function(e) GetAssayData(obj, assay = "RNA", slot = layer))
if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0)
  stop(sprintf("layer '%s' 为空", layer))
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "CsparseMatrix")

# --- 诊断：所选层数值范围（确认 counts=整数 / data≈log归一化 max<~15） ---
cat(sprintf("export layer = %s\n", layer))
xv <- mat@x[seq_len(min(3000, length(mat@x)))]
cat(sprintf("%s: min=%g max=%g integer=%s\n", layer, min(xv), max(xv), all(xv == round(xv))))

Matrix::writeMM(mat, file.path(outdir, "counts.mtx"))   # 文件名固定，build_h5ad.py 直接读入 X
writeLines(rownames(mat), file.path(outdir, "genes.txt"))
writeLines(colnames(mat), file.path(outdir, "barcodes.txt"))
write.csv(obj@meta.data, file.path(outdir, "meta.csv"))
cat(sprintf("exported %d genes x %d cells (layer=%s) -> %s\n", nrow(mat), ncol(mat), layer, outdir))
