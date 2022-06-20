
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

scRNA = readRDS(file = "keyobject/fbsubid.rds")

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- FindClusters(scRNA, resolution = 0.1)

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)

table(scRNA$fbsub)

table(Idents(scRNA))

##作图

dir.create("harmony/cluster2")

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cluster2/UMAP.png", plot = plot2, width = 5.5, height = 5)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cluster2/tsne.png", plot = plot2, width = 5.5, height = 5)

plot2 = DimPlot(scRNA, reduction = "umap", group.by = "fbtype")

ggsave("harmony/cluster2/UMAPfbt.png", plot = plot2, width = 6, height = 5)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cluster2/markers.rds")

all.markers = readRDS(file = "harmony/cluster2/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cluster2/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cluster2/top50_diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "harmony/cluster2/top20mk.csv", row.names = F)

features = c(top50$gene)

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

ggsave("harmony/cluster2/top50mk.png", p, width = 10, height = 15)

submk = subset(top50, cluster == 2)

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk,
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cluster2/submkvio.png", width = 25, height = 4.5)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

fbsub = data.frame(ClusterID = 0:(n - 1), fbsub = "unkown")  #构建数据框

table(Idents(scRNA))

## 判断亚群ID属于那类细胞

fbsub[fbsub$ClusterID %in% c(1, 7), 2] = "FB_1"

fbsub[fbsub$ClusterID %in% c(2), 2] = "FB_2"

fbsub[fbsub$ClusterID %in% c(3), 2] = "FB_3"

fbsub[fbsub$ClusterID %in% c(4), 2] = "FB_4"

fbsub[fbsub$ClusterID %in% c(5), 2] = "FB_5"

fbsub[fbsub$ClusterID %in% c(6, 0), 2] = "FB_6"

## 重新赋值

scRNA@meta.data$fbsub = "NA"

for(i in 1:nrow(fbsub)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == fbsub$ClusterID[i]), "fbsub"] <- fbsub$fbsub[i]}

table(scRNA@meta.data$fbsub)

table(scRNA$seurat_clusters)

dir.create("harmony/cellsub2022")

Idents(scRNA) = scRNA@meta.data$fbsub

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cellsub2022/UMAP.png", plot = plot2, width = 8, height = 7)

plot2 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3)

ggsave("harmony/cellsub2022/tsne.png", plot = plot2, width = 6, height = 5)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cellsub2022/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cellsub2022/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cellsub2022/top50_diff_genes_wilcox.csv", row.names = F)

top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

features = c(top50$gene)

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

ggsave("harmony/cellsub2022/top50mk.png", p, width = 6, height = 10)

pal <- paletteer_d("ggsci::category20_d3")[c(3, 6, 4, 5, 2, 1)]

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

submk = subset(top20, cluster == "FB_1")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio1.png", width = 8, height = 3.5)

submk = subset(top20, cluster == "FB_2")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio2.png", width = 8, height = 3.5)

submk = subset(top20, cluster == "FB_3")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio3.png", width = 8, height = 3.5)

submk = subset(top20, cluster == "FB_4")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio4.png", width = 8, height = 3.5)

submk = subset(top20, cluster == "FB_5")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio5.png", width = 8, height = 3.5)

submk = subset(top20, cluster == "FB_6")

mk = submk$gene

mk

p1 = VlnPlot(scRNA, features = mk, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "harmony/cellsub2022/submkvio6.png", width = 8, height = 3.5)

table(scRNA$fbsub)

n = length(unique(scRNA@meta.data$seurat_clusters)) #获取亚群个数

fbsub = data.frame(ClusterID = 0:(n - 1), fbsub = "unkown")  #构建数据框

## 判断亚群ID属于那类细胞

fbsub[fbsub$ClusterID %in% c(1, 7), 2] = "DPT_FB"

fbsub[fbsub$ClusterID %in% c(2), 2] = "BAMBI_FB"

fbsub[fbsub$ClusterID %in% c(3), 2] = "CLDN1_FB"

fbsub[fbsub$ClusterID %in% c(4), 2] = "RGS5_FB"

fbsub[fbsub$ClusterID %in% c(5), 2] = "CXCL14_FB"

fbsub[fbsub$ClusterID %in% c(6, 0), 2] = "COL11A1_FB"

## 重新赋值

scRNA@meta.data$fbsub = "NA"

for(i in 1:nrow(fbsub)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == fbsub$ClusterID[i]), "fbsub"] <- fbsub$fbsub[i]}

table(scRNA@meta.data$fbsub)

table(scRNA$seurat_clusters)

dir.create("harmony/cafsub2022")

Idents(scRNA) = scRNA@meta.data$fbsub

plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3)

ggsave("harmony/cafsub2022/UMAP.png", plot = plot2, width = 6, height = 5)

plan()

diff.wilcox = FindAllMarkers(scRNA)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val < 0.05)

saveRDS(all.markers, file = "harmony/cafsub2022/markers.rds")

all.markers = readRDS(file = "harmony/cafsub2022/markers.rds")

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "harmony/cafsub2022/top10_diff_genes_wilcox.csv", row.names = F)

top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(top50, "harmony/cafsub2022/top50_diff_genes_wilcox.csv", row.names = F)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "harmony/cafsub2022/top20mk.csv", row.names = F)

features = c(top50$gene)

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

ggsave("harmony/cafsub2022/top50mk.png", p, width = 10, height = 15)

table(scRNA$fbsub)

table(Idents(scRNA))

saveRDS(scRNA, file = "keyobject/fbsubid.rds")
