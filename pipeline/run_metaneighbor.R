#!/usr/bin/env Rscript
# run_metaneighbor.R — canonical MetaNeighbor (Bioconductor) 跨队列复现性 AUROC。
# 用法: Rscript run_metaneighbor.R [scratch_dir]   (默认 ./work/tmp/mn)
# 读 <scratch_dir>/{expr.mtx, genes.txt, meta.csv} → 写 <scratch_dir>/auroc.csv
#   expr.mtx : genes × cells (log 归一化即可，MetaNeighbor 内部按秩处理)
#   meta.csv : 列 cohort, indep_cl（行序与 expr 列序一致）
# 安装: R -e 'BiocManager::install("MetaNeighbor")'  (没 BiocManager 先 install.packages)
suppressMessages({
  library(MetaNeighbor); library(SummarizedExperiment); library(Matrix)
})
a <- commandArgs(trailingOnly = TRUE)
d <- if (length(a) >= 1) a[1] else "./work/tmp/mn"
expr <- as(Matrix::readMM(file.path(d, "expr.mtx")), "CsparseMatrix")   # genes × cells
genes <- readLines(file.path(d, "genes.txt"))
rownames(expr) <- genes
meta <- read.csv(file.path(d, "meta.csv"), colClasses = "character")
stopifnot(ncol(expr) == nrow(meta))
colnames(expr) <- paste0("c", seq_len(ncol(expr)))

se <- SummarizedExperiment(assays = list(expr = expr),
                           colData = DataFrame(study_id = meta$cohort,
                                               cell_type = meta$indep_cl))
cat(sprintf("MetaNeighbor 输入: %d genes × %d cells; %d 队列; %d 个(队列|簇)节点\n",
            nrow(expr), ncol(expr), length(unique(meta$cohort)),
            length(unique(paste(meta$cohort, meta$indep_cl)))))

vg <- variableGenes(dat = se, exp_labels = meta$cohort)
cat(sprintf("variableGenes: %d\n", length(vg)))

auroc <- MetaNeighborUS(var_genes = vg, dat = se,
                        study_id = meta$cohort, cell_type = meta$indep_cl,
                        fast_version = TRUE)
write.csv(auroc, file.path(d, "auroc.csv"))
cat(sprintf("AUROC 矩阵 %d×%d → %s/auroc.csv\n", nrow(auroc), ncol(auroc), d))
cat(sprintf("AUROC 分布: min=%.3f median=%.3f max=%.3f\n",
            min(auroc, na.rm = TRUE), median(auroc, na.rm = TRUE), max(auroc, na.rm = TRUE)))
