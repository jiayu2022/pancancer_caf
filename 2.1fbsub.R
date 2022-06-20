
getwd()

library(Seurat)

library(tidyverse)

library(patchwork)

library(dplyr)

library("ggsci")

library(paletteer)

library(pheatmap) # 加载包

library(MySeuratWrappers)

library(Scillus)

rm(list = ls())

#dir.create("figure")

#dir.create("figure/fbsub")

scRNA = readRDS(file = "keyobject/fbsubid.rds")

table(scRNA$fbsub)

table(Idents(scRNA))

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "fbsub", label = F, repel = T, cols = pal) +
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = "")) +
        theme(legend.key.height = unit(0.8, "cm"))

ggsave("figure/fbsub/umap1c20.png", plot = plot1, width = 5.4, height = 4)

ggsave("figure/fbsub/umap1c20.pdf", plot = plot1, width = 5.4, height = 4)

pal <- paletteer_d("ggsci::category20_d3")[c(20:1)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "cell_orig_ident", label = F, repel = T, cols = pal) +
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = "")) +
        theme(legend.key.height = unit(0.5, "cm"))

ggsave("figure/fbsub/umap2.png", plot = plot1, width = 5.7, height = 4)

ggsave("figure/fbsub/umap2.pdf", plot = plot1, width = 5.7, height = 4)

pal <- paletteer_d("ggsci::category20_d3")[c(2, 1)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "fbtype", label = F, repel = T, cols = pal) +
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = "")) +
        theme(legend.key.height = unit(0.5, "cm"))

ggsave("figure/fbsub/umapfbt.png", plot = plot1, width = 5.4, height = 4)

ggsave("figure/fbsub/umapfbt.pdf", plot = plot1, width = 5.4, height = 4)

pal <- paletteer_d("ggsci::category20_d3")[c(1, 2, 5)]

plot1 = DimPlot(scRNA, reduction = "umap", group.by = "Phase", label = F, repel = T, cols = pal) +
        ggtitle(NULL) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + #加边框
        theme(axis.line.x = element_line(colour = "black", size = 0.1), axis.line.y = element_line(colour = "black", size = 0.1)) +
        guides(fill = guide_legend(title = "")) +
        theme(legend.key.height = unit(1, "cm"))

ggsave("figure/fbsub/umapcc.png", plot = plot1, width = 5.1, height = 4)

ggsave("figure/fbsub/umapcc.pdf", plot = plot1, width = 5.1, height = 4)

scRNA$seurat_clusters = scRNA$fbsub

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "sampletype", pal_setup = pal) +
            theme_bw() + coord_flip() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "", x = "") +
            guides(fill = guide_legend(title = "Cluster"))

ggsave("figure/fbsub/cellfreq.png", plot1, width = 6, height = 4)

ggsave("figure/fbsub/cellfreq.pdf", plot1, width = 6, height = 4)

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "sampletype", pal_setup = pal) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "", x = "") +
            guides(fill = guide_legend(title = "Cluster"))

ggsave("figure/fbsub/cellfreqxs.png", plot1, width = 9, height = 6)

ggsave("figure/fbsub/cellfreqxs.pdf", plot1, width = 9, height = 6)

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "fbtype", pal_setup = pal) +
            theme_bw() + coord_flip() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "", x = "") +
            guides(fill = guide_legend(title = "Cluster"))

ggsave("figure/fbsub/fbfreq.png", plot1, width = 8, height = 2)

ggsave("figure/fbsub/fbfreq.pdf", plot1, width = 8, height = 2)

plot1 = plot_stat(scRNA, plot_type = "prop_fill", group_by = "fbtype", pal_setup = pal) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(y = "", x = "") +
            guides(fill = guide_legend(title = "Cluster"))

ggsave("figure/fbsub/fbfreqxs.png", plot1, width = 3, height = 6)

ggsave("figure/fbsub/fbfreqxs.pdf", plot1, width = 3, height = 6)

plan()

table(Idents(scRNA))

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

top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

write.csv(top100, "figure/marker/top100deg.csv", row.names = F)

top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

write.csv(top15, "figure/marker/top15deg.csv", row.names = F)

top15 = read.csv(file = "figure/marker/top15deg.csv")

features = c(top15$gene)

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = features,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

#dir.create("harmony/fbsub")

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

#mypal <- pal_d3("category10", alpha = 0.8)(10)

#mypal = c("#1F77B4CC", "#FF7F0ECC", "#2CA02CCC", "#D62728CC", "#9467BDCC",
#          "#8C564BCC", "#E377C2CC", "#7F7F7FCC", "#BCBD22CC", "#17BECFCC")

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "",
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

ggsave("figure/marker/top15mk.png", p, width = 3.5, height = 8)

ggsave("figure/marker/top15mk.pdf", p, width = 3.5, height = 8)

table(scRNA$fbsub)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

markers <- c("PDGFRA",
             "FAP", "PDPN",
             "LUM", "COL6A3", "THY1", "DCN",
             "COL3A1", "COL1A1", "COL1A2")

p1 = VlnPlot(scRNA, features = markers, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker/fbmkvio.png", width = 5, height = 3.5)

ggsave(p1, file = "figure/marker/fbmkvio.pdf", width = 5, height = 3.5)

markergene = c("PDGFRA",
               "FAP", "PDPN",
               "LUM", "COL6A3", "THY1", "DCN",
               "COL3A1", "COL1A1")

p2 <- FeaturePlot(scRNA, features = markergene, reduction = "umap", label = F, label.size = 3, ncol = 3,
                  cols = c("lightgrey", "red"))

ggsave("figure/marker/fbmkumap.png", p2, width = 10, height = 9)

ggsave("figure/marker/fbmkumap.pdf", p2, width = 10, height = 9)

#dir.create("harmony/fbsub")

annotation_col = data.frame(Cluster = factor(c("BAMBI+", "CLDN1+", "COL11A1+", "CXCL14+", "DPT+", "RGS5+" )))

rownames(annotation_col)

colnames(mfdf)

rownames(annotation_col) <- colnames(mfdf)

head(annotation_col)

#mypal <- pal_d3("category10", alpha = 0.8)(10)

#mypal = c("#1F77B4CC", "#FF7F0ECC", "#2CA02CCC", "#D62728CC", "#9467BDCC",
#          "#8C564BCC", "#E377C2CC", "#7F7F7FCC", "#BCBD22CC", "#17BECFCC")

ann_colors = list(Cluster = c("DPT+" = "#2CA02CCC", "CLDN1+" = "#8C564BCC", "CXCL14+" = "#D62728CC",
                              "BAMBI+" = "#9467BDCC", "COL11A1+" = "#FF7F0ECC", "RGS5+" = "#1F77B4CC"))

head(ann_colors)

ECM = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2",
        "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL10A1",
        "COL11A1", "COL12A1", "COL14A1",
        "FN1", "POSTN", "MMP2", "MMP9", "MMP14", "TIMP1",
        "TIMP2", "TIMP3", "LOX")

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = ECM,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "ECM-associated",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker/ecm.png", p, width = 3, height = 8)

ggsave("figure/marker/ecm.pdf", p, width = 3, height = 8)

IM = c("IL1B", "IL6", "IL8", "IL10",
       "CXCL1", "CXCL5", "CXCL9", "CXCL10", "CXCL12", "CXCL14",
       "CCL2", "CCL5", "CCL7",
       "CSF1",
       "CD274", "PDCD1LG2", "CD276", "NT5E",
       "DPP4", "JAM2", "TNFSF4", "CHI3L1", "PTGES", "TDO2")

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = IM,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "Immune Response-associated",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker/im.png", p, width = 4, height = 8.5)

ggsave("figure/marker/im.pdf", p, width = 4, height = 8.5)

GF = c("FGF1", "FGF2", "FGF7",
       "IGF1", "IGF2", "CTGF", "HGF",
       "VEGFA", "VEGFB", "VEGFC",
       "PDGFA", "PDGFC")

mk <- AverageExpression(scRNA,
                        group.by = "fbsub",
                        features = GF,
                        return.seurat = FALSE)

mfdf = as.data.frame(mk[["RNA"]])

p <- pheatmap(mfdf,
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              annotation_legend = FALSE,
              annotation_names_col = FALSE,
              color = colorRampPalette(c("#0080FF", "white", "#FF8000"))(20),
              main = "Growth Factor",
              scale = "row",
              angle_col = 270, # 设置显示角度
              cellwidth = 20, cellheight = 20, # 设置热图方块宽度和高度
              border = "white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 10, # 分别设置横向和纵向字体大小
              fontsize_col = 10.5)

ggsave("figure/marker/gf.png", p, width = 3, height = 5)

ggsave("figure/marker/gf.pdf", p, width = 3, height = 5)

pal <- paletteer_d("ggsci::category20_d3")[c(5, 6, 2, 4, 3, 1)]

ECM = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2",
        "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL10A1",
        "COL11A1", "COL12A1", "COL14A1",
        "FN1", "POSTN", "MMP2", "MMP9", "MMP14", "TIMP1",
        "TIMP2", "TIMP3", "LOX")

p1 = VlnPlot(scRNA, features = ECM, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker/ecmvio.png", width = 10, height = 3.5)

ggsave(p1, file = "figure/marker/ecmvio.pdf", width = 10, height = 3.5)

IM = c("IL1B", "IL6", "IL8", "IL10",
       "CXCL1", "CXCL5", "CXCL9", "CXCL10", "CXCL12", "CXCL14",
       "CCL2", "CCL5", "CCL7",
       "CSF1",
       "CD274", "PDCD1LG2", "CD276", "NT5E",
       "DPP4", "JAM2", "TNFSF4", "CHI3L1", "PTGES", "TDO2")

p1 = VlnPlot(scRNA, features = IM, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker/imvio.png", width = 10, height = 3.5)

ggsave(p1, file = "figure/marker/imvio.pdf", width = 10, height = 3.5)

GF = c("FGF1", "FGF2", "FGF7",
       "IGF1", "IGF2", "CTGF", "HGF",
       "VEGFA", "VEGFB", "VEGFC",
       "PDGFA", "PDGFC")

p1 = VlnPlot(scRNA, features = GF, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank())  #不显示坐标刻度

ggsave(p1, file = "figure/marker/gfvio.png", width = 6, height = 3.5)

ggsave(p1, file = "figure/marker/gfvio.pdf", width = 6, height = 3.5)

IM = c("IL1B", "IL6", "IL8", "IL10",
       "CXCL1", "CXCL5", "CXCL9", "CXCL10", "CXCL12", "CXCL14",
       "CCL2", "CCL5", "CCL7",
       "CSF1",
       "CD274", "PDCD1LG2", "CD276", "NT5E",
       "DPP4", "JAM2", "TNFSF4", "CHI3L1", "PTGES", "TDO2")

p1 = VlnPlot(scRNA, features = IM, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker/imvioa.png", width = 10, height = 3.5)

ggsave(p1, file = "figure/marker/imvioa.pdf", width = 10, height = 3.5)

GF = c("FGF1", "FGF2", "FGF7",
       "IGF1", "IGF2", "CTGF", "HGF",
       "VEGFA", "VEGFB", "VEGFC",
       "PDGFA", "PDGFC")

p1 = VlnPlot(scRNA, features = GF, cols = pal, group.by = "fbsub",
             stacked = T, pt.size = 0,
             direction = "horizontal", #水平作图
             x.lab = "", y.lab = "") + #横纵轴不标记任何东西
        theme(axis.text.x = element_text(angle = 90, size = rel(0.35))) +
        theme(axis.ticks.x = element_line(size = 0.35))

ggsave(p1, file = "figure/marker/gfvioa.png", width = 6, height = 3.5)

ggsave(p1, file = "figure/marker/gfvioa.pdf", width = 6, height = 3.5)
