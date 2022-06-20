
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library(harmony)

rm(list = ls())

library(future)

plan()

options(future.globals.maxSize = 32 * 1024^3)

plan("multicore", workers = 56)

scRNA = readRDS(file = "keyobject/scRNAcafnf.rds")

#dir.create("harmony")

dir.create("harmony/qc")

dir.create("harmony/cluster")

scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

p2 <- VlnPlot(scRNA, features = "percent.rb", group.by = "cell_orig_ident", pt.size = 0.01) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("harmony/qc/rbqcb.png", p2, width = 8, height = 5)

scRNA <- subset(scRNA, subset = percent.rb < 40)

p2 <- VlnPlot(scRNA, features = "percent.rb", group.by = "cell_orig_ident", pt.size = 0.01) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("harmony/qc/rbqca.png", p2, width = 8, height = 5)

scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData()

s.genes <- cc.genes$s.genes

g2m.genes <- cc.genes$g2m.genes

scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

head(scRNA[[]])

scRNA <- RunPCA(scRNA, features = c(s.genes, g2m.genes))

plot1 = DimPlot(scRNA)

ggsave("harmony/cluster/pcacc.png", plot = plot1, width = 8, height = 6.5)

scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))

scRNA <- RunPCA(scRNA, features = c(s.genes, g2m.genes))

plot1 = DimPlot(scRNA)

ggsave("harmony/cluster/pcaccsd.png", plot = plot1, width = 8, height = 6.5)

scRNA <- RunPCA(scRNA, verbose = FALSE)

#dir.create("tmp")

saveRDS(scRNA, file = "tmp/fbqca.rds")

scRNA <- RunHarmony(scRNA, group.by.vars = "batch")

#降维聚类

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- FindClusters(scRNA, resolution = 0.1)

##作图

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cluster/UMAP.png", plot = plot2, width = 8, height = 7)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cell_orig_ident")

ggsave("harmony/cluster/UMAP1.png", plot = plot1, width = 8, height = 6.5)

plot2 = DimPlot(scRNA, reduction = "umap", group.by = "fbtype")

ggsave("harmony/cluster/UMAPfbt.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cluster/tsne.png", plot = plot2, width = 8, height = 7)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by = "cell_orig_ident")

ggsave("harmony/cluster/tsne1.png", plot = plot1, width = 8, height = 6.5)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by = "fbtype")

ggsave("harmony/cluster/tsne2.png", plot = plot1, width = 8, height = 6.5)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "Phase")

ggsave("harmony/cluster/umapcc.png", plot = plot1, width = 8, height = 6.5)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by = "Phase")

ggsave("harmony/cluster/tsnecc.png", plot = plot1, width = 8, height = 6.5)

select_genes <- c("EPCAM", "PTPRC", "PECAM1", "CSPG4", "COL1A1", "ACTA2")

#featureplot展示

p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "tsne", label = F, label.size = 3, ncol = 3,
                  cols = c("lightgrey", "red"))

ggsave("harmony/cluster/celltype_marker.png", p2, width = 10, height = 6)

select_genes <- c("ACTA2", "COL1A1", "FAP", "THY1", "PDPN", "S100A4")

#featureplot展示
p2 <- FeaturePlot(scRNA, features = select_genes, reduction = "tsne", label = F, label.size = 3, ncol = 3,
                  cols = c("lightgrey", "red"))

ggsave("harmony/cluster/celltype_marker5f.png", p2, width = 10, height = 6)

plan()

table(Idents(scRNA))

#Idents(scRNA) = scRNA$seurat_clusters

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cluster/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cluster/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cluster/top50_diff_genes_wilcox.csv", row.names = F)

top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

features = c(top100$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "CAFsub Markers",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20
              )

ggsave("harmony/cluster/top100mk.png", p, width = 10, height = 15)

table(scRNA$seurat_clusters)

table(scRNA@meta.data$fbtype)

Cells.sub <- subset(scRNA@meta.data, seurat_clusters != 7 & seurat_clusters != 8 &
                                     seurat_clusters != 9 & seurat_clusters != 10)

scRNAsub <- subset(scRNA, cells = row.names(Cells.sub))

table(scRNAsub$seurat_clusters)

saveRDS(scRNAsub, file = "tmp/scRNAsubcaf.rds")

scRNA = scRNAsub

scRNA = readRDS(file = "tmp/scRNAsubcaf.rds")

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)

dir.create("harmony/cellsub")

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cellsub/UMAP.png", plot = plot2, width = 8, height = 7)

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cell_orig_ident")

ggsave("harmony/cellsub/UMAP1.png", plot = plot1, width = 8, height = 6.5)

plot2 = DimPlot(scRNA, reduction = "umap", group.by = "fbtype")

ggsave("harmony/cellsub/UMAPfbt.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cellsub/tsne.png", plot = plot2, width = 8, height = 7)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by = "cell_orig_ident")

ggsave("harmony/cellsub/tsne1.png", plot = plot1, width = 8, height = 6.5)

plot1 = DimPlot(scRNA, reduction = "tsne", group.by = "fbtype")

ggsave("harmony/cellsub/tsne2.png", plot = plot1, width = 8, height = 6.5)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cellsub/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cellsub/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cellsub/top50_diff_genes_wilcox.csv", row.names = F)

top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

features = c(top100$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "CAFsub Markers",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20
              )

ggsave("harmony/cellsub/top100mk.png", p, width = 10, height = 15)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

fbsub = data.frame(ClusterID = 0:(n - 1), fbsub = "unkown")  #构建数据框

## 判断亚群ID属于那类细胞

fbsub[fbsub$ClusterID %in% c(1), 2] = "FB_1"

fbsub[fbsub$ClusterID %in% c(2), 2] = "FB_2"

fbsub[fbsub$ClusterID %in% c(3), 2] = "FB_3"

fbsub[fbsub$ClusterID %in% c(4, 0), 2] = "FB_4"

fbsub[fbsub$ClusterID %in% c(5), 2] = "FB_5"

fbsub[fbsub$ClusterID %in% c(6), 2] = "FB_6"

## 重新赋值

scRNA@meta.data$fbsub = "NA"

for(i in 1:nrow(fbsub)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == fbsub$ClusterID[i]), "fbsub"] <- fbsub$fbsub[i]}

table(scRNA@meta.data$fbsub)

table(scRNA$seurat_clusters)

#dir.create("harmony/cellsub2")

Idents(scRNA) = scRNA@meta.data$fbsub

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cellsub2/UMAP.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cellsub2/tsne.png", plot = plot2, width = 8, height = 7)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cellsub2/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cellsub2/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cellsub2/top50_diff_genes_wilcox.csv", row.names = F)

top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

features = c(top100$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "CAFsub Markers",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20
              )

ggsave("harmony/cellsub2/top100mk.png", p, width = 10, height = 15)

table(scRNA$fbsub)

plot2 = DimPlot(scRNA, reduction = "umap", group.by = "Phase")

ggsave("harmony/cellsub2/UMAPcc.png", plot = plot2, width = 8, height = 7)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

fbsub = data.frame(ClusterID = 0:(n - 1), fbsub = "unkown")  #构建数据框

## 判断亚群ID属于那类细胞

fbsub[fbsub$ClusterID %in% c(1), 2] = "RGS5_FB"

fbsub[fbsub$ClusterID %in% c(2), 2] = "CFD_FB"

fbsub[fbsub$ClusterID %in% c(3), 2] = "HP_FB"

fbsub[fbsub$ClusterID %in% c(4, 0), 2] = "POSTN_FB"

fbsub[fbsub$ClusterID %in% c(5), 2] = "CXCL14_FB"

fbsub[fbsub$ClusterID %in% c(6), 2] = "CCL2_FB"

## 重新赋值

scRNA@meta.data$fbsub = "NA"

for(i in 1:nrow(fbsub)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == fbsub$ClusterID[i]), "fbsub"] <- fbsub$fbsub[i]}

table(scRNA@meta.data$fbsub)

table(scRNA$seurat_clusters)

dir.create("harmony/cafsub")

Idents(scRNA) = scRNA@meta.data$fbsub

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cafsub/UMAP.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cafsub/tsne.png", plot = plot2, width = 8, height = 7)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cafsub/markers.rds")

all.markers = readRDS(file = "harmony/cafsub/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cafsub/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cafsub/top50_diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "harmony/cafsub/top20mk.csv", row.names = F)

features = c(top20$gene)

library(pheatmap) # 加载包

mk <- AverageExpression(
  scRNA,
  features = features,
  return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              main = "CAFsub Markers",
              scale = "row",
              angle_col = 270, # 设置显示角度
              #cellwidth = 20,cellheight = , # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 1, # 分别设置横向和纵向字体大小
              fontsize_col = 20
              )

ggsave("harmony/cafsub/top100mk.png", p, width = 10, height = 15)

table(scRNA$fbsub)

plot2 = DimPlot(scRNA, reduction = "umap", group.by = "Phase")

ggsave("harmony/cafsub/UMAPcc.png", plot = plot2, width = 8, height = 7)

table(Idents(scRNA))

saveRDS(scRNA, file = "keyobject/fbsubid.rds")

scRNA = readRDS(file = "keyobject/fbsubid.rds")

scRNA <- FindVariableFeatures(scRNA) %>% ScaleData() %>% RunPCA(verbose = FALSE)

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)

plot2 = DimPlot(scRNA, reduction = "umap", label = F, label.size = 3)

ggsave("harmony/cafsub/UMAPa.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = F, label.size = 3)

ggsave("harmony/cafsub/tsnea.png", plot = plot2, width = 8, height = 7)

sampletype <- c("2020_CRC_Adj" = "CRC_Adj",
                "2020_LC_Adj" = "LC_Adj",
                "2020_OVC_Adj" = "OVC_Adj",
                "2020_ICC_Adj" = "ICC_Adj",
                "2020_PDAC_Adj" = "PDAC_Adj",
                "2021_OVC_Adj" = "OVC_Adj",
                "2019_HCC_ICC" = "HCC_ICC",
                "2020_BC" = "BC",
                "2020_CRC" = "CRC",
                "2020_ICC" = "ICC",
                "2020_LC" = "LC",
                "2020_OVC" = "OVC",
                "2020_PDAC" = "PDAC",
                "2021_ICC" = "ICC",
                "2021_OVC" = "OVC",
                "2021_PRAD" = "PRAD")

scRNA$sampletype = unname(sampletype[scRNA$cell_orig_ident])

table(scRNA$sampletype)

saveRDS(scRNA, file = "keyobject/fbsubid.rds")
