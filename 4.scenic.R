
getwd()

#dir.create("figure/scenic")

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library("ggsci")

library(paletteer)

library(pheatmap) # 加载包

#library(SCENIC)

library(AUCell)

rm(list = ls())

library(future)

plan()

options(future.globals.maxSize = 32 * 1024^3)

plan("multisession", workers = 56)

scRNA = readRDS(file = "keyobject/fbsubid.rds")

table(Idents(scRNA))

###### 导出pyscenic需要使用的表达矩阵

mat.all <- as.matrix(scRNA@assays$RNA@data)

# 过滤表达量过低的基因，减少后续分析的计算量

mat.filter <- mat.all[rowSums(mat.all) > 200, ]

write.csv(t(mat.filter), file = "figure/scenic/caf_all.csv")

str(scRNA)

#pyscenic

regulonAUC <- read.csv("figure/scenic/auc_mtx.csv", row.names = 1, check.names = F)

regulonAUC <- t(regulonAUC)

dir.create("figure/scenic/process")

saveRDS(regulonAUC, file = "figure/scenic/process/regulonAUC.rds")

regulonAUC = readRDS(file = "figure/scenic/process/regulonAUC.rds")

bin.T <- AUCell_exploreThresholds(regulonAUC,
                                  smallestPopPercent = 0.25,
                                  assignCells = TRUE,
                                  plotHist = FALSE,
                                  verbose = FALSE)

regulonBin <- lapply(rownames(regulonAUC), function(reg){
  as.numeric(colnames(regulonAUC) %in% bin.T[[reg]][["assignment"]])
})

regulonBin <- do.call("rbind", regulonBin)

dimnames(regulonBin) <- list(rownames(regulonAUC), colnames(regulonAUC))

saveRDS(regulonBin, "figure/scenic/process/regulonBin.rds")

regulonBin = readRDS(file = "figure/scenic/process/regulonBin.rds")

metadata = scRNA@meta.data

scAUC = CreateSeuratObject(regulonAUC,
                           project = "AUC",
                           meta.data = metadata)

scBin = CreateSeuratObject(regulonBin,
                           project = "Bin",
                           meta.data = metadata)

saveRDS(scAUC, "figure/scenic/process/scAUC.rds")

saveRDS(scBin, "figure/scenic/process/scBin.rds")

scBin = readRDS(file = "figure/scenic/process/scBin.rds")

scRNA = scBin

Idents(scRNA) <- scRNA@meta.data$fbsub

levels(scRNA) #查看是否已改名

#默认wilcox方法

diff.wilcox = FindAllMarkers(scRNA, logfc.threshold = 0.2, only.pos = T)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "figure/scenic/process/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(all.markers, "figure/scenic/process/diff_genes_wilcox.csv", row.names = F)

write.csv(top10, "figure/scenic/process/top10_diff_genes_wilcox.csv", row.names = F)

top6 = all.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)

write.csv(top6, "figure/scenic/process/top6_diff_genes_wilcox.csv", row.names = F)

top10 = read.csv(file = "figure/scenic/top10a.csv")

features = c(top10$gene)

table(scRNA$fbsub)

mk <- AverageExpression(scRNA,
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
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "SCENIC",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 5, # 分别设置横向和纵向字体大小
              fontsize_col = 8)

#dir.create("figure/marker")

ggsave("figure/scenic/top10tf.png", p, width = 3, height = 6)

ggsave("figure/scenic/top10tf.pdf", p, width = 3, height = 6)

