
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(monocle)

library("ggsci")

library(paletteer)

library(future)

plan()

options(future.globals.maxSize = 32 * 1024^3)

plan("multicore", workers = 56)

###

rm(list = ls())

scRNA = readRDS(file = "keyobject/fbsubid.rds")

table(Idents(scRNA))

sessionInfo()

data <- GetAssayData(scRNA, assay = "RNA", slot = "counts")

cell_metadata <- scRNA@meta.data

colnames(cell_metadata)

table(cell_metadata$sampletype)

cell_meta_sub = select(cell_metadata, sampletype, fbtype, Phase, fbsub)

colnames(cell_meta_sub)

gene_annotation <- data.frame(gene_short_name = rownames(data))

rownames(gene_annotation) <- rownames(data)

pd <- new("AnnotatedDataFrame", data = cell_meta_sub)

fd <- new("AnnotatedDataFrame", data = gene_annotation)

cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
                      phenoData = pd, featureData = fd,
                      expressionFamily = negbinomial.size(),
                      lowerDetectionLimit = 1)

class(cds)

cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds)

dim(cds)

dir.create("figure/monocle")

saveRDS(cds, file = "figure/monocle/cds_raw.rds")

cds <- detectGenes(cds, min_expr = 0.1)

disp_table <- dispersionTable(cds)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

cds <- reduceDimension(cds,
                       max_components = 2,
                       num_dim = 30,
                       reduction_method = 'tSNE', verbose = T, check_duplicates = FALSE)

cds <- clusterCells(cds, num_clusters = 6)

p1 = plot_cell_clusters(cds, 1, 2,
                        color = "fbsub")

ggsave("figure/monocle/cluster1.png", p1, width = 5, height = 6)

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~fbsub",
                                      cores = 56)

saveRDS(cds, file = "figure/monocle/cds1.rds")

saveRDS(diff_test_res, file = "figure/monocle/diff_test_res.rds")

cds = readRDS(file = "figure/monocle/cds1.rds")

diff_test_res = readRDS(file = "figure/monocle/diff_test_res.rds")

sig_genes <- subset(diff_test_res, qval < 0.1)

dim(sig_genes)

head(sig_genes[, c("gene_short_name", "pval", "qval")])

cg = head(sig_genes$gene_short_name)

p = plot_genes_jitter(cds[cg, ],
                      grouping = "fbsub",
                      ncol = 2)

ggsave("figure/monocle/plotg.png", p, width = 6, height = 6)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) ; length(ordering_genes)

cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2,
                       reduction_method = "DDRTree")

# 细胞排序

cds <- orderCells(cds)

# 可视化

p = plot_cell_trajectory(cds, color_by = "fbsub")

ggsave("figure/monocle/ddrtree.png", p, width = 5, height = 5)

ggsave("figure/monocle/ddrtree.pdf", p, width = 5, height = 5)

# 展现marker基因

p = plot_genes_in_pseudotime(cds[cg, ],
                             color_by = "fbsub")

ggsave("figure/monocle/ddrgene.png", p, width = 5, height = 5)

ggsave("figure/monocle/ddrgene.pdf", p, width = 5, height = 5)
