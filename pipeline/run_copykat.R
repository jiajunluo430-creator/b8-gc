#!/usr/bin/env Rscript
# run_copykat.R — 对单个样本的上皮 counts 跑 CopyKAT，输出 aneuploid/diploid 预测。
# 用法: Rscript run_copykat.R <in_prefix> <out_dir> <sample> [n_cores]
#   in_prefix: <scratch_dir>/in/<sample>  (读 .mtx/.genes/.barcodes)
#   依赖: copykat, Matrix    安装: R -e 'remotes::install_github("navinlabcode/copykat")'
#
# 输出 <out_dir>/<sample>_copykat_prediction.txt （CopyKAT 默认在工作目录生成，脚本会搬运）
suppressMessages({library(copykat); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
prefix <- a[1]; outdir <- a[2]; sample <- a[3]
ncores <- if (length(a) >= 4) as.integer(a[4]) else 8
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

mat <- as(Matrix::readMM(paste0(prefix, ".mtx")), "CsparseMatrix")     # genes × cells
genes <- readLines(paste0(prefix, ".genes"))
cells <- readLines(paste0(prefix, ".barcodes"))
rownames(mat) <- genes; colnames(mat) <- cells
cat(sprintf("[%s] %d genes × %d cells\n", sample, nrow(mat), ncol(mat)))

wd <- getwd(); setwd(outdir)
res <- tryCatch(
  copykat(rawmat = as.matrix(mat), id.type = "S", ngene.chr = 5,
          win.size = 25, KS.cut = 0.1, sam.name = sample,
          distance = "euclidean", norm.cell.names = "", output.seg = FALSE,
          plot.genes = FALSE, genome = "hg20", n.cores = ncores),
  error = function(e) {cat(sprintf("[%s] CopyKAT 失败: %s\n", sample, conditionMessage(e))); NULL})
setwd(wd)

if (!is.null(res)) {
  pred <- data.frame(res$prediction)
  write.table(pred, file.path(outdir, paste0(sample, "_copykat_prediction.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] 完成: %s\n", sample,
              paste(names(table(pred$copykat.pred)), table(pred$copykat.pred),
                    sep = "=", collapse = " ")))
}
