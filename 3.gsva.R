
getwd()

library(Seurat)

library(tidyverse)

library(GSVA)

library(msigdbr)

library("ggsci")

library(paletteer)

library(patchwork)

library(pheatmap) # 加载包

#dir.create("figure/gsva")

#dir.create("figure/gsva/tmp")

rm(list = ls())

scRNA = readRDS(file = "keyobject/fbsubid.rds")

table(Idents(scRNA))

library(future)

plan()

options(future.globals.maxSize = 32 * 1024^3)

plan("multicore", workers = 56)

#hallmarker

genesets <- msigdbr(species = "Homo sapiens", category = "H") %>% select("gs_name","gene_symbol") %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

str(genesets)

table(Idents(scRNA))

expr <- as.matrix(scRNA@assays$RNA@counts)

gsvago = gsva(expr, genesets, method = "ssgsea", parallel.sz = 56)

saveRDS(gsvago, file = "figure/gsva/tmp/gsvah.rds")

gsvago = readRDS(file = "figure/gsva/tmp/gsvah.rds")

str(gsvago)

gsvago[1:5, 1:10]

str(scRNA@meta.data)

meta <- scRNA@meta.data[, c("seurat_clusters", "fbsub")]

meta[1:5, 1:2]

scg = CreateSeuratObject(gsvago, project = "Hallmark", meta.data = meta)

table(Idents(scg))

Idents(scg) = scg@meta.data$fbsub

table(Idents(scg))

saveRDS(scg, file = "figure/gsva/tmp/sch.rds")

scg = readRDS(file = "figure/gsva/tmp/sch.rds")

library(patchwork)

str(scg)

str(scg@assays$RNA@data@Dimnames[1])

hm = scg@assays$RNA@data@Dimnames[1]

str(hm)

features = c(hm[[1]])

features

mk <- AverageExpression(scg,
                        group.by = "fbsub",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = TRUE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = TRUE,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/gsva/htop15mka.png", p, width = 4.5, height = 6)

ggsave("figure/gsva/htop15mka.pdf", p, width = 4.5, height = 6)

table(scRNA$fbsub)

#GO-BP

m_df = msigdbr(species = "Homo sapiens")

m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)

genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% select("gs_name","gene_symbol") %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

str(genesets)

table(Idents(scRNA))

scRNAsub = subset(scRNA, downsample = 1000)

table(Idents(scRNAsub))

saveRDS(scRNAsub, file = "figure/gsva/tmp/scRNAsub1000.rds")

plan()

expr <- as.matrix(scRNAsub@assays$RNA@counts)

gsvago = gsva(expr, genesets, method = "ssgsea", parallel.sz = 56)

saveRDS(gsvago, file = "figure/gsva/tmp/gsvago.rds")

gsvago = readRDS(file = "figure/gsva/tmp/gsvago.rds")

str(gsvago)

gsvago[1:5, 1:10]

str(scRNAsub@meta.data)

meta <- scRNAsub@meta.data[, c("seurat_clusters", "fbsub")]

meta[1:5, 1:2]

scg = CreateSeuratObject(gsvago, project = "GO_BP", meta.data = meta)

table(Idents(scg))

Idents(scg) = scg@meta.data$fbsub

table(Idents(scg))

saveRDS(scg, file = "figure/gsva/tmp/scg.rds")

#scg = readRDS(file = "figure/gsva/tmp/scg.rds")

plan()

table(Idents(scg))

diff.wilcox = FindAllMarkers(scg, logfc.threshold = 0.02)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

#dir.create("figure/gsva/go")

saveRDS(all.markers, file = "figure/gsva/go/markers.rds")

#all.markers = readRDS(file = "figure/gsva/go/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "figure/gsva/go/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "figure/gsva/go/top50_diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "figure/gsva/go/top20mk.csv", row.names = F)

#没必要

goset = read.csv(file = "figure/gsva/go/top15deg.csv")

features = c(goset$gene)

features

mk <- AverageExpression(scg,
                        group.by = "fbsub",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = TRUE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = TRUE,
              treeheight_col = 20, # 分别设置横、纵向聚类树高
              treeheight_row = 20,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

ggsave("figure/gsva/go/gotop15mka.png", p, width = 4.5, height = 6)

ggsave("figure/gsva/go/gotop15mka.pdf", p, width = 4.5, height = 6)

table(scRNA$fbsub)

